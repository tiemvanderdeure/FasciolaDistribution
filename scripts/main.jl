using FasciolaDistribution
using Rasters, RasterDataSources, Proj, NCDatasets, StatsBase, 
    Random, Turing, Tables, Statistics, Dates, Rasters.Lookups
import CSV
const FD = FasciolaDistribution # any function called with FD. is from this package

# this is where the data will be stored
ENV["FASCIOLA_DATA_PATH"] = joinpath(RasterDataSources.rasterpath(), "FasciolaDistribution")

isdir(FD.datapath()) || mkdir(FD.datapath())

global chunks = (X(256), Y(256))
global force = false
global ncwritekw = (; chunks, force = true, deflatelevel = 2, shuffle = false)

# Set some global variables
# Region of interest
using NaturalEarth
countries = naturalearth("admin_0_countries", 10)
europe = countries.geometry[countries.CONTINENT .== "Europe" .&& countries.ADM0_A3 .!= "RUS"]
africa = countries.geometry[countries.CONTINENT .== "Africa"]
roi = [europe; africa]
ext = Extent(X=(-25, 53), Y = (-35, 73)) # roughly Africa + Europe

##############################
# Thermal performance curves #
##############################
# Fit temperature-dependent performance of life history traits related to fasciola transmission based on experimental data
# Most objects in this section share a nested NamedTuple structure. The outer NamedTuple has 
# keys `hepatica` and `gigantica` and the inner has keys corresponding to traits (see FasciolaDistribution.traits_symbol)

# Load in experimental life history data. 
life_history_data = FD.load_life_history_data()

# Define models. This binds data to the models without sampling from them
life_history_models = FD.define_life_history_models(life_history_data)

# Sample from the models. 
life_history_chains = map(life_history_models) do models
    map(models) do m
        chain = sample(Xoshiro(0), m, NUTS(), 1000, progress=true)
    end
end

# Each model returns a Vector of callable structs as the generated quantites
# Each struct contains a function and the estimated parameters at that sample. 
# They can be called with just temperature as input and return the value of the life history trait
gqs = map(life_history_models, life_history_chains) do models, chains
    map(models, chains) do m, c
        vec(generated_quantities(m, c))
    end
end

# For convenience later on, convert all the namedtuple of vectors to a vector of namedtuples

# Use Tables.rowtable to convert to a vector of NamedTuples
gqs_rt = map(Tables.rowtable, gqs)

# `gqs_rt.hepatica`is a vector with as element type a NamedTuple of callable structs, 
# corresponding to each life history trait

###################################
# Thermal suitability Projections #
###################################
# Use the fits from above to project the thermal suitability of fasciola transmission

# run this and save, to save lots of time.
if force || !isfile(joinpath(datapath(), "temperature_current.nc"))
    temperature = FD.get_temperature_data(ext; ssps, gcms, dates)

    map((:current, :future)) do t
        path = joinpath(datapath(), "temperature_$(t).nc")
        write(path, FD.writeable_dims(temperature[t]); ncwritekw...)
    end
end

# read in temperature rasters, then fix dimensions
temperature = NamedTuple(t => Raster(joinpath(datapath(), "temperature_$(t).nc")) for t in (:current, :future))
temperature = map(temperature) do ras
    setcrs(
        set(ras, dims((gcms, ssps, dates), dims(ras))...),
        EPSG(4326)
    )
end

### convert temperature into suitability
# map over current and future
temperature_suitability = map(temperature) do t_avg
    # get unique temperatures
    unique_ts = sort(collect(keys(StatsBase.countmap(t_avg)))) # for some reason much faster than Base.unique

    # convert temperature to a Raster of indices
    # each index refers to a temperature in unique_ts
    indices = searchsortedfirst.(Ref(sort(unique_ts)), t_avg);
    missingidx = findfirst(Base.Fix2(===, missingval(t_avg)), unique_ts)

    # map over species
    transmission_rasters = map(gqs_rt) do g
        # run the life cycle model for unique temperatures
        transmissionstrength = life_cycle_model.(unique_ts', g)
        # normalize so all values are between 0 and 1
        transmissionstrength ./= maximum(view(transmissionstrength, :, Not(missingidx)), dims = 2)

        mean_transmissionstrength = dropdims(mean(transmissionstrength; dims = 1); dims = 1)
        mean_transmissionstrength[missingidx] = NaN # this is the missing value

        # do this with a view to save memory
        transmission_raster = rebuild(indices; 
            data = view(mean_transmissionstrength, parent(indices)), 
            missingval = NaN
        ) 
    end
end

# takes ~2 minutes
FD.write_rasters(joinpath(datapath(), "monthly_temperature_suitability"), temperature_suitability; ncwritekw...)

temperature = nothing
GC.gc()

### Write posterior estimates for optimal transmission temperature to a text file
open("images/optimal_temp_posterior.txt", "w") do io
    map((:hepatica, :gigantica)) do s
        ts = 10:0.01:40
        hep_suitability = FD.life_cycle_model.(collect(ts)', gqs_rt[s])
        optimal_temps = getindex.(Ref(ts), last.(findmax.(eachrow(hep_suitability))))
        println(io, "optimal temp $s $(mean(optimal_temps)), 95% CI $(quantile(optimal_temps, [0.025, 0.975]))")
    end
end

########################################
# Hydrological suitability Projections #
########################################

# Read in discharge data and save to disk
if force || !isfile(joinpath(datapath(), "discharge_current.nc"))
    discharge = get_discharge_data(ext; gcms, ssps, ghms, dates)
    map((:current, :future)) do t
        path = joinpath(datapath(), "discharge_$(t).nc")
        write(path, writeable_dims(discharge[t]); ncwritekw...)
    end
end

# Read data from disk
discharge = map((:current, :future)) do t
    ras = Raster(joinpath(datapath(), "discharge_$(t).nc"); raw = true)
    setcrs(
        set(ras, dims((ghms, ssps, dates), dims(ras))...),
        EPSG(4326)
    )    
end |> NamedTuple{(:current, :future)}

# Normalize so units are mm/day instead of m3/s
m3_s_to_mm_day = (1000*60*60*24) ./ cellarea(first(discharge))
discharge_mm_day = map(discharge) do x
    x .* m3_s_to_mm_day
end;

# Disaggregate to the same resolution as the temperature data
discharge_res = map(discharge_mm_day) do d
    # use lazy (view) to save memory
    da = disaggregate(d, (X(12), Y(12)); lazy = true)
    # make sure dims are identical to avoid nonsense later
    set(da, dims(temperature_suitability.current.hepatica))
end

# For now, assume 50% hydrological suitability at runoff of 1 mm/day
discharge_to_suitability(x) = x / (one(x) + x)
hydrological_suitability = map(x -> (FD.@lazyd discharge_to_suitability.(x)), discharge_res)

# takes ~5 minutess
FD.write_rasters(joinpath(datapath(), "monthly_hydro_suitability"), hydrological_suitability; ncwritekw...)

##################################
# Combined transmission strength #
##################################
# Get monthly transmission suitability by multiplying temperature and hydrological suitability for each month
monthly_transmission_suitability = map(temperature_suitability, hydrological_suitability) do ts, h
    map(ts) do t
        @d FD.geomean.(h, t) strict=false
    end
end

# takes ~20 minutes
FD.write_rasters(joinpath(datapath(), "monthly_parasite_suitability"), monthly_transmission_suitability; ncwritekw...)

# Takes ~30 minutes
for file in filter(startswith("monthly"), readdir(datapath()))
    ras = Raster(joinpath(datapath(), file); lazy = true)
    newfile = replace(file, "monthly_" => "")
    open(ras) do x
        annual = FD.f_and_dropdims(mean, x; dims = :month)
        summarized = FD.f_and_dropdims(mean, annual; dims = Rasters.commondims(ras, (:ghm, :gcm)))
        @info "writing $newfile"
        write(joinpath(datapath(), newfile), annual; ncwritekw...)
        write(joinpath(datapath(), "mean_" * newfile), summarized; ncwritekw...)
    end
end

#####################
# Vector suitabilty #
#####################
# Load Maxnet library
using Maxnet, StatisticalMeasures
import SpeciesDistributionModels as SDM

# First read in bioclim data
if force || !isfile(joinpath(datapath(), "bio_current.nc"))
    predictors = Tuple(i for i in (1:19) if !(i in (8,9,18,19)))
    bio = get_bioclim(predictors, ext; ssps, gcms, dates) # replace missing  in here?
    for t in (:current, :future)
        path = joinpath(datapath(), "bio_$(t).nc")
        write(path, writeable_dims(bio[t]); ncwritekw...)
    end
end

# Read data from disk
bio = map((:current, :future)) do t
    ras = RasterStack(joinpath(datapath(), "bio_$(t).nc"); raw = true)
    ras = set(ras, dims((ghms, ssps, dates), dims(ras))...)
    setcrs(ras, EPSG(4326))
end |> NamedTuple{(:current, :future)}

predictors = bio

# Read in snail occurrence data
# This includes some basic cleaning like selecting on year and removing duplicates
occurrences, samplingbg = FD.get_snail_occurrences()

# occurrences have galba split into europe and africa. All radix is Africa
rois = (galba_eu = europe, galba_af = africa, radix = roi)

# Further processing - extract bioclimatic data, thin to one sample per grid cell,
# sample random background points, increase the weight of African Galba.
occ_bgs = map(keys(occurrences)) do K
    o, sbg = map((occurrences, samplingbg)) do x
        # extract from raster
        e = extract(predictors.current, x[K], skipmissing = true, index = true, geometry = true)
        # find unique indices (= grid cells)
        unique_indices = unique(i -> e[i].index, eachindex(e))
        # select only unique grid cells
        e = e[unique_indices]
        # split into geometry and bioclimatic data - geometry is used just for plotting
        geo = getfield.(e, :geometry)
        bio_e = Base.structdiff.(e, NamedTuple{(:index,:geometry)})
        return (; geo, bio = bio_e)
    end
    # sample sampling background points up to the number of presences points
    n_ocs = length(o.geo)
    n_bg = length(sbg.geo)
    if n_bg > n_ocs
        indices = sample(eachindex(sbg.geo), n_ocs; replace = false)
        sbg = map(x -> x[indices], sbg)
    end
    
    # draw random background points from within the region
    b = mask(predictors.current; with = roi)
    n_unif = n_ocs
    unif = Rasters.sample(
        Xoshiro(0), b, n_unif; weights = cellarea(b), geometry = (X,Y), skipmissing = true
        # Xoshiro(0) is for reproducibility
    )
    # split into geography and predictors
    ubg = (
        geo = getfield.(unif, :geometry),
        bio = Base.structdiff.(unif, NamedTuple{(:geometry,)})
    )
    # combine sampling bg points and uniformly sampled bg points
    bg = map(vcat, sbg, ubg)

    return (o, bg)
end |> NamedTuple{keys(occurrences)}

# Finally put african and european galba together
bgs = (galba = map(vcat, occ_bgs.galba_eu[2], occ_bgs.galba_af[2]), radix = occ_bgs.radix[2])
ocs = (galba = map(vcat, occ_bgs.galba_eu[1], occ_bgs.galba_af[1]), radix = occ_bgs.radix[1])

open("images/number_occurrences.txt","w") do io
    map(keys(ocs)) do K
        println(io, "total $K: $(length(ocs[K].geo))")
    end
    println(io, "Galba Africa: $(length(occ_bgs.galba_af[1].geo))")
end

#=
# Predictor selection
biomat = Tables.matrix(vcat(ocs.galba.bio, bgs.galba.bio))
cmat = cor(biomat)
biokeys = keys(predictors.current)
cmat_nt = NamedTuple{biokeys}(NamedTuple{biokeys}.(eachrow(cmat)))
to_include = (:bio1, :bio7, :bio12, :bio17)
candidates = filter(x -> all(y -> abs(x[y]) < 0.7, to_include), cmat_nt)
=#
models = (; maxnet = MaxnetBinaryClassifier(features = "lqp", regularization_multiplier = 3))

evs = map(ocs, bgs) do o, bg
    cvdata = SDM.sdmdata(o.bio, bg.bio; 
        predictors = (:bio1, :bio7, :bio12, :bio17), resampler = SDM.CV(rng = Xoshiro(0), nfolds = 5))
    m = SDM.sdm(cvdata, models)
    measures = (auc = SDM.auc, tss = BalancedAccuracy(adjusted = true))
    ev = SDM.evaluate(m; measures)
    mach_evs = SDM.machine_evaluations(ev)
    mapreduce(vcat, keys(mach_evs)) do k
        mapreduce(vcat, keys(mach_evs[k])) do k2
            (; dataset = k, measure = k2, scores = mean(mach_evs[k][k2]))
        end
    end
end

for k in keys(evs)
    CSV.write(joinpath("images", "evaluation_$(k).csv"), evs[k])
end

# fit models on all data
models = map(ocs, bgs) do o, bg
    data = SDM.sdmdata(o.bio, bg.bio; 
        predictors = (:bio1, :bio7, :bio12, :bio17))
    m = SDM.sdm(data, models)
end

# ## Predict suitability for both species and for current and future climate
# Map over models (species)
host_suitability = map(predictors) do b
    # Map over bioclimatic data (current and future)
    map(models) do m
        dropdims(SDM.predict(m, b); dims = :Band)    
    end
end;

FD.write_rasters(joinpath(datapath(), "host_suitability"), host_suitability; ncwritekw...)

#################
# Save to files #
#################

# Get annual transmission suitability taking the mean over months
# The actual workhorse - should take ~ 10 minutes
for t in (:current, :future)
    @show t
    for (f, snail) in ((:hepatica, :galba), (:gigantica, :radix))
        @show f
        transmission_suitability = Raster(joinpath(datapath(), "transmission_suitability_$(t)_$(f).nc"); lazy = true)
        host_suitability = Raster(joinpath(datapath(), "host_suitability_$(t)_$(snail).nc"); lazy = true)
        # this is the final projection
        fasciola_risk = @d FD.geomean.(transmission_suitability, host_suitability) strict=false
        fasciola_risk = set(fasciola_risk, dims(host_suitability))
        open(fasciola_risk) do x
            write(
                joinpath(datapath(), "fasciola_risk_$(t)_$(f).nc"),
                x;
                ncwritekw...
            )
        end
    end
end


##### Calculate means over GCMs/GHMs for plotting
# Read everything in lazily
lazypreds =
    FD.read_predictions((gcms, ghms, dates); lazy = true, raw = true)

# get the filenames
filenames = joinpath.(datapath(), "mean_" .* (
    "hydrological_suitability",
    "host_suitability",
    "transmission_suitability",
    "temperature_suitability",
    "fasciola_risk"
    )
)

# summarize and save each file
map(lazypreds, filenames) do r, file
    FD.write_rasters(file, r; ncwritekw...) do x
        # exclude this GCM as it has very high ECS
        x2 = hasdim(x, :gcm) ? view(x, gcm = Not(At(UKESM1_0_LL))) : x
        FD.f_and_dropdims(mean, x2; dims = Rasters.commondims(x2, (:gcm, :ghm)))
    end
end

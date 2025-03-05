using FasciolaDistribution
using Rasters, RasterDataSources, Proj, NCDatasets, StatsBase, 
    Random, Turing, Tables, Statistics, Dates, Rasters.Lookups
const FD = FasciolaDistribution # any function called with FD. is from this package

# this is where the data will be stored
ENV["FASCIOLA_DATA_PATH"] = joinpath(RasterDataSources.rasterpath(), "FasciolaDistribution")

isdir(datapath()) || mkdir(datapath())

global chunks = (X(256), Y(256))
global force = false
global ncwritekw = (; chunks, force = true, deflatelevel = 2, shuffle = false)

# Set some global variables
# Region of interest
# TODO: What is our region of interest??
using NaturalEarth
countries = naturalearth("admin_0_countries", 10)
europe = countries.geometry[countries.CONTINENT .== "Europe" .&& countries.ADM0_A3 .!= "RUS"]
africa = countries.geometry[countries.CONTINENT .== "Africa"]
roi = [europe; africa]
ext = Extent(X=(-25.0001, 53.0001), Y = (-35.0001, 73.0001))



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

temperature = nothing
GC.gc()

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
    # make sure dims are identical to aviod nonsense later
    #set(da, dims(temperature_suitability.current.hepatica))
end

# For now, assume 50% hydrological suitability at runoff of 1 mm/day
discharge_to_suitability(x) = x / (one(x) + x)
hydrological_suitability = map(x -> (FD.@lazyd discharge_to_suitability.(x)), discharge_res)

##################################
# Combined transmission strength #
##################################

# Get monthly transmission suitability by multiplying temperature and hydrological suitability for each month
monthly_transmission_suitability = map(temperature_suitability, hydrological_suitability) do ts, h
    map(ts) do t
        FD.@lazyd t .* h
    end
end

#####################
# Vector suitabilty #
#####################
# Load Maxnet library
using Maxnet

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

waterdist = FD.water_distance(bio.current)

# combine bioclimatic and geographic predictors
predictors = map(bio) do b
    merge(b, (; waterdist))
end

predictors = bio#map(b -> b[(:bio1, :bio7, :bio12, :bio15)], bio)

# Read in snail occurrence data
# This includes some basic cleaning like selecting on year and removing duplicates
occurrences, samplingbg = FD.get_snail_occurrences()

# occurrences have galba split into europe and africa. All radix is Africa
weights = (galba_eu = 1, galba_af = 3, radix = 1)
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
    # increase the weight of African Galba by repeating the entries
    o = map(x -> repeat(x, weights[K]), o)
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

# ## Fit the models
models = map(ocs, bgs) do o, bg
    y = vcat(trues(length(o.bio)), falses(length(bg.bio)))
    x = vcat(o.bio, bg.bio)
    maxnet(y, x; features = "lqp", regularization_multiplier = 2)
end

# ## Predict suitability for both species and for current and future climate
host_suitability = map(models) do m
    mypredict(m, predictors.current)
end;

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])
plot!(ax1, host_suitability.galba, colorrange = (0,1), colormap = Reverse(:Spectral))
plot!(ax2, host_suitability.radix, colorrange = (0,1), colormap = Reverse(:Spectral))
fig

# Map over models (species)
host_suitability = map(predictors) do b
    # Map over bioclimatic data (current and future)
    map(models) do m
        mypredict(m, b)
    end
end;

bio = nothing; GC.gc()

for t in (:current, :future)
    for snail in (:galba, :radix)
        write(
        joinpath(datapath(), "host_suitability_$(t)_$(snail).nc"),
        FD.writeable_dims(host_suitability[t][snail]);
        ncwritekw...
    )
    end
end

hydro, host, trans, temp, risk =
    FD.read_predictions((gcms, ghms, dates); lazy = true)
host_suitability = host


#################
# Save to files #
#################

# Get annual transmission suitability taking the mean over months
# The actual workhorse - should take ~ 10 minutes
for t in (:current, :future)
    @show t
    for (f, snail) in ((:hepatica, :galba), (:gigantica, :radix))
        @show f
        transmission_suitability = f_and_dropdims(mean, monthly_transmission_suitability[t][f]; dims = :month)
        write(
            joinpath(datapath(), "transmission_suitability_$(t)_$(f).nc"),
            FD.writeable_dims(transmission_suitability);
            ncwritekw...
        )
        
        write(
            joinpath(datapath(), "temperature_suitability_$(t)_$f.nc"),
            FD.writeable_dims(f_and_dropdims(mean, temperature_suitability[t][f]; dims = :month));
            ncwritekw...
        )
            
        write(
            joinpath(datapath(), "host_suitability_$(t)_$(snail).nc"),
            FD.writeable_dims(host_suitability[t][snail]);
            ncwritekw...
        )
        # this is the final projection
        fasciola_risk = @d sqrt.(transmission_suitability .* host_suitability[t][snail]) strict=false
        write(
            joinpath(datapath(), "fasciola_risk_$(t)_$(f).nc"),
            FD.writeable_dims(fasciola_risk);
            ncwritekw...
        )
    end

    write(
        joinpath(datapath(), "hydrological_suitability_$(t).nc"),
        FD.writeable_dims(dropdims(mean(hydrological_suitability[t]; dims = :month); dims = :month));
        ncwritekw...
    )
end

### Calculate the mean for all of these objects so they fit in memory
f_and_dropdims(f, x; dims) = dropdims(f(x; dims); dims)

function mean_or_loop(x, filename)
    if x isa Raster
        println("writing $filename")
        open(x) do x
            ras = f_and_dropdims(mean, x; dims = otherdims(x, (X,Y,:ssp,Ti)))
            write(
                filename * ".nc", FD.writeable_dims(ras);
                ncwritekw...
            )
        end
    elseif x isa NamedTuple
        for K in keys(x)
            mean_or_loop(x[K], filename * "_" * string(K))
        end
    end
end

filenames = joinpath.(datapath(), "mean_" .* (
    "hydrological_suitability",
    "host_suitability",
    "transmission_suitability",
    "temperature_suitability",
    "fasciola_risk"
    )
)

hydro, host, trans, temp, risk =
    FD.read_predictions((gcms, ghms, dates); lazy = true, raw = true)

mean_or_loop(host, joinpath(datapath(), "mean_host_suitability"))    

# lazily calculate mean for all of these
for (x, filename) in zip((hydro, host, trans, temp, risk), filenames)
    mean_or_loop(x, filename)    
end



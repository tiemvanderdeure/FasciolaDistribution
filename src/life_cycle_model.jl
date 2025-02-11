function life_cycle_model(
    x,
    hatchingtime,
    hatchingsuccess,
    egg_death_rate,
    infectionefficiency,
    prepatentperiod,
    snail_death_rate,
    cercarialrelease,
)
    ht = hatchingtime(x)
    isfinite(ht) || return 0.0
 
    hs = hatchingsuccess(x)
    iszero(hs) && return 0.0

    edr = egg_death_rate(x)
    isfinite(edr) || return 0.0

    ie = infectionefficiency(x)
    iszero(ie) && return 0.0

    pp = prepatentperiod(x)
    isfinite(pp) || return 0.0

    sdr = snail_death_rate(x)
    isfinite(sdr) || return 0.0

    cr = cercarialrelease(x)
    iszero(cr) && return 0.0

    hs * exp(-edr * ht) * ie * exp(-sdr * pp) * cr
end
function life_cycle_model(x, fs::NamedTuple{trait_keys})
    life_cycle_model(x, fs.hatching_time, fs.hatching_success, x -> 0.01,
        fs.infection_efficiency, fs.prepatent_period, x -> 0.01, fs.cercarial_release)
end
life_cycle_model(::Missing, fs::NamedTuple{trait_keys}) = missing

function find_max_life_cycle(gqs)
    to_optimize(x) = -life_cycle_model(x[1], gqs)
    -Optim.minimum(optimize(to_optimize, [20.0]))[1]
end
#=
include("fits.jl")
using Rasters, RasterDataSources, NaturalEarth, ArchGDAL

countries = naturalearth("admin_0_countries", 10)
africa = countries.geometry[countries.CONTINENT .== "Africa"]

tas = RasterSeries(CHELSA{Climate}, :tas; month = 1:12, lazy = true)
tas_af = read(crop(tas; to = africa))
tas_af = cat(tas_af..., dims = Ti(1:12))
tas_af_m = mask(tas_af; with = africa)
md = Rasters.metadata(tas_af_m)
tas_normalized = tas_af_m .* md["scale"] .+ md["offset"]

bm = boolmask(tas_normalized)
tas_sm = collect(skipmissing(tas_normalized))

unique(tas_af_m)

plot(@view(tas_af_m[Ti(1)]))

@profview res = life_cycle_model.(
    collect(0.0:0.001:30.0)', 
    gqs_hepatica.hatchingtime, 
    Ref(hatchingsuccess), 
    Ref(egg_death_rate), 
    gqs_hepatica.infectionefficiency, 
    gqs_hepatica.prepatentperiod, 
    Ref(snail_death_rate), 
    gqs_hepatica.cercarialrelease
)

x = 19.0
@time [life_cycle_model(
    x, 
    gqs_hepatica
) for x in 1.0:10.0]

@time gqs_hepatica.hatchingtime[i].(1.0:100.0)

@time myfun.(10.0:20.0, Ref(gqs_gigantica.hatchingtime[1]))

function myfun(x, fix)
    fix(x)
end


function life_cycle_model(x, gqs, i = 1)
    life_cycle_model(
    x, 
    gqs.hatchingtime[i], 
    x -> 0.5, 
    x -> 0.01, 
    gqs.infectionefficiency[i], 
    gqs.prepatentperiod[i], 
    x -> 0.01, 
    gqs.cercarialrelease[i]
)
end

function life_cycle_model(
    t, 
    tmin_hatchingtime, gdd_mu_hatchingtime, 
    rmax_infectionefficiency, infection_efficiency_tmin, infection_efficiency_tmax, 
    tmin_prepatentperiod, gdd_mu_prepatentperiod,
    rmax_cercarialrelease, Topt_cercarialrelease, breadth_cercarialrelease
)
    if t < tmin_hatchingtime || t < tmin_prepatentperiod || t < infection_efficiency_tmin || t > infection_efficiency_tmax
        0.0
    else
        egg_death_rate = 0.01
        egg_hatching_time = gdd_time(tmin_hatchingtime, gdd_mu_hatchingtime, t)
        infection_efficiency = quadratic(rmax_infectionefficiency, infection_efficiency_tmin, infection_efficiency_tmax, t)
        prepatent_period = gdd_time(tmin_prepatentperiod, gdd_mu_prepatentperiod, t)
        snail_death_rate = 0.01
        metacercarial_release = exp.(gaussian(rmax_cercarialrelease, Topt_cercarialrelease, breadth_cercarialrelease, t))

        exp(-egg_death_rate * egg_hatching_time) * infection_efficiency * exp(-snail_death_rate * prepatent_period) * metacercarial_release
    end
end

function life_cycle_model_gig(
    t, 
    tmin_hatchingtime, gdd_hatchingtime, 
    tmin_hatchingsuccess, tmax_hatchingsuccess, a_hatchingsuccess,
    rmax_infectionefficiency, infection_efficiency_tmin, infection_efficiency_tmax, 
    tmin_prepatentperiod, gdd_mu_prepatentperiod,
    rmax_cercarialrelease, Topt_cercarialrelease, breadth_cercarialrelease
)
    if t < tmin_hatchingtime || t < tmin_prepatentperiod || t < infection_efficiency_tmin || t > infection_efficiency_tmax ||
        t < tmin_hatchingsuccess
        0.0
    else
        egg_death_rate = 0.01
        egg_hatching_time = gdd_time(tmin_hatchingtime, gdd_hatchingtime, t)
        hatching_success = logistic_tmin(tmin_hatchingsuccess, tmax_hatchingsuccess, a_hatchingsuccess, t)
        infection_efficiency = quadratic(rmax_infectionefficiency, infection_efficiency_tmin, infection_efficiency_tmax, t)
        prepatent_period = gdd_time(tmin_prepatentperiod, gdd_mu_prepatentperiod, t)
        snail_death_rate = 0.01
        metacercarial_release = exp.(gaussian(rmax_cercarialrelease, Topt_cercarialrelease, breadth_cercarialrelease, t))

        hatching_success * exp(-egg_death_rate * egg_hatching_time) * infection_efficiency * exp(-snail_death_rate * prepatent_period) * metacercarial_release
    end
end

life_cycle_model(::Missing, args...) = missing

tmin_hatchingtime = chains_hepatica.hatchingtime[:Tmin]
gdd_hatchingtime = chains_hepatica.hatchingtime[:gdd_mu]

tmin_prepatentperiod = chains_hepatica.prepatentperiod[:Tmin]
gdd_prepatentperiod = chains_hepatica.prepatentperiod[:gdd_mu]

width_infectionefficiency = chains_hepatica.infectionefficiency[:width]
rmax_infectionefficiency = chains_hepatica.infectionefficiency[:rmax]
Topt_infectinonefficiency = chains_hepatica.infectionefficiency[:Topt]
infection_efficiency_tmin = Topt_infectinonefficiency .- width_infectionefficiency
infection_efficiency_tmax = Topt_infectinonefficiency .+ width_infectionefficiency

tmin_prepatentperiod = chains_hepatica.prepatentperiod[:Tmin]
gdd_prepatentperiod = chains_hepatica.prepatentperiod[:gdd_mu]

rmax_cercarialrelease = chains_hepatica.cercarialrelease[:rmax_mu]
Topt_cercarialrelease = chains_hepatica.cercarialrelease[:Topt]
breadth_cercarialrelease = chains_hepatica.cercarialrelease[:breadth]

temps = 0.0:0.01:40.0

@time @fastmath metacercariae_per_egg = life_cycle_model.(
    temps', 
    tmin_hatchingtime[1:1000], gdd_hatchingtime[1:1000], 
    rmax_infectionefficiency, infection_efficiency_tmin, infection_efficiency_tmax, 
    tmin_prepatentperiod[1:1000], gdd_prepatentperiod[1:1000], 
    rmax_cercarialrelease, Topt_cercarialrelease, breadth_cercarialrelease
);

rel_metacercariae_per_egg = metacercariae_per_egg ./ maximum(metacercariae_per_egg; dims = 2)
replace!(rel_metacercariae_per_egg, NaN => 0.0)

med_vec_cap = median(rel_metacercariae_per_egg; dims = 1)
quant_vec_cap = reduce(hcat, quantile.(eachcol(rel_metacercariae_per_egg), Ref([0.025,0.975])))

maxtemp = temps[getindex.(findmax(metacercariae_per_egg; dims = 2)[2], 2)]
quantile(maxtemp, [0.025,0.975])

fig, ax, l = lines(temps, med_vec_cap[1,:])
for i in 1:2
    lines!(ax, temps, quant_vec_cap[i,:], linestyle = :dash)
end
fig


### gigantica
tmin_hatchingtime = chains_gigantica.hatchingtime[:Tmin]
gdd_hatchingtime = chains_gigantica.hatchingtime[:gdd_mu]

tmin_hatchingsuccess = chains_gigantica.hatchingsuccess[:Tmin]
tmax_hatchingsuccess = chains_gigantica.hatchingsuccess[:Tmax]
a_hatchingsuccess = chains_gigantica.hatchingsuccess[:a]

tmin_prepatentperiod = chains_gigantica.prepatentperiod[:Tmin]
gdd_prepatentperiod = chains_gigantica.prepatentperiod[:gdd_mu]

width_infectionefficiency = chains_gigantica.infectionefficiency[:width]
rmax_infectionefficiency = chains_gigantica.infectionefficiency[:rmax]
Topt_infectinonefficiency = chains_gigantica.infectionefficiency[:Topt]
infection_efficiency_tmin = Topt_infectinonefficiency .- width_infectionefficiency
infection_efficiency_tmax = Topt_infectinonefficiency .+ width_infectionefficiency

tmin_prepatentperiod = chains_gigantica.prepatentperiod[:Tmin]
gdd_prepatentperiod = chains_gigantica.prepatentperiod[:gdd_mu]

rmax_cercarialrelease = chains_gigantica.cercarialrelease[:rmax_mu]
Topt_cercarialrelease = chains_gigantica.cercarialrelease[:Topt]
breadth_cercarialrelease = chains_gigantica.cercarialrelease[:breadth]

@time @fastmath metacercariae_per_egg = life_cycle_model_gig.(
    temps', 
    tmin_hatchingtime, gdd_hatchingtime, 
    tmin_hatchingsuccess, tmax_hatchingsuccess, a_hatchingsuccess,
    rmax_infectionefficiency, infection_efficiency_tmin, infection_efficiency_tmax, 
    tmin_prepatentperiod, gdd_prepatentperiod, 
    rmax_cercarialrelease, Topt_cercarialrelease, breadth_cercarialrelease
);

rel_metacercariae_per_egg = metacercariae_per_egg ./ maximum(metacercariae_per_egg; dims = 2)
replace!(rel_metacercariae_per_egg, NaN => 0.0)

med_vec_cap = median(rel_metacercariae_per_egg; dims = 1)
quant_vec_cap = reduce(hcat, quantile.(eachcol(rel_metacercariae_per_egg), Ref([0.025,0.975])))

maxtemp = temps[getindex.(findmax(metacercariae_per_egg; dims = 2)[2], 2)]
quantile(maxtemp, [0.025,0.975])
median(maxtemp)

fig, ax, l = lines(temps, med_vec_cap[1,:])
for i in 1:2
    lines!(ax, temps, quant_vec_cap[i,:], linestyle = :dash)
end
fig





dict = Dict(unique_temps .=> metacercariae_per_egg[1,:])
get(x, dict) = dict[x]
get(::Missing, dict) = missing

using BenchmarkTools
vals =tas_normalized[100_000_000:110_000_000]
vals_sm = collect(skipmissing(vals))
@btime get.($vals, Ref($dict))

@profview vc = get.(@view(tas_normalized[Ti(1)]), Ref(dict))

unique_temps = unique(tas_normalized)
unique_temps = collect(skipmissing(unique_temps))

@time metacercariae_per_egg = life_cycle_model.(
    tas_normalized[Ti(1)], 
    tmin_hatchingtime[1], gdd_hatchingtime[1], 
    rmax_infectionefficiency[1], infection_efficiency_tmin[1], infection_efficiency_tmax[1], 
    tmin_prepatentperiod[1], gdd_prepatentperiod[1], 
    rmax_cercarialrelease[1], Topt_cercarialrelease[1], breadth_cercarialrelease[1]
)

plot(metacercariae_per_egg)

plot(tas_normalized[Ti(1)])
=#
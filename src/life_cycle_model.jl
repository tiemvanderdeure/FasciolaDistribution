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
function life_cycle_model(
    x, # temperature 
    fs::NamedTuple{trait_keys}; # namedtuple with structs with trait functions
    egg_death_rate = 0.01, # default egg death rate
    snail_death_rate = 0.01 # default snail death rate
)
    life_cycle_model(x, fs.hatching_time, fs.hatching_success, x -> egg_death_rate,
        fs.infection_efficiency, fs.prepatent_period, x -> snail_death_rate, fs.cercarial_release)
end
life_cycle_model(::Missing, fs::NamedTuple{trait_keys}; kw...) = missing

function find_max_life_cycle(gqs)
    to_optimize(x) = -life_cycle_model(x[1], gqs)
    -Optim.minimum(optimize(to_optimize, [20.0]))[1]
end
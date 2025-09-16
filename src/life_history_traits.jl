### Life history traits
# Functions to read in the life history data and to define the Bayesian
# models to fit thermal performance curves to these

## Loading data
isvalidunit(unit) = in(unit, ("days", "N", "cercariae/snail")) || Base.contains("radioactivity")(unit)

function load_life_history_data()
    data = CSV.read(joinpath(@__DIR__, "..", "data", "life_history.csv"), DataFrame)
    data_incl = filter(row -> row.include && isvalidunit(row.unit), data)
    disallowmissing!(data_incl, [:temperature, :observation])
    data_gig = filter(row -> !ismissing(row.fasciola) && row.fasciola == "gigantica", data_incl)
    data_hep = filter(row -> !ismissing(row.fasciola) && row.fasciola == "hepatica", data_incl)
    data_natalensis = filter(row -> ismissing(row.fasciola) && row.snail in ["radix natalensis", "Lymnaea cailliaudi"], data_incl)
    data_truncatula = filter(row -> ismissing(row.fasciola) && row.snail == "Galba truncatula", data_incl)
    
    data_hep_galba = filter(row -> ismissing(row.snail) || row.snail == "Galba truncatula", data_hep)
    # for now, include radix auricularia as well. Without it, there are huge gaps in the data (e.g. prepatent period at low temperatures)
    data_gig_radix = filter(
        row -> ismissing(row.snail) || row.snail == "radix natalensis" || 
            (Base.contains("auricularia")(row.snail) && row.trait == "prepatent period"), 
        data_gig
    )

    map((hepatica = data_hep_galba, gigantica = data_gig_radix)) do data
        map(traits_str) do trait
            filter(row -> row.trait == trait, data)
        end |> NamedTuple{trait_keys}
    end
end

function life_history_citations()
    lifehistorybib = BibParser.parse_file(joinpath(@__DIR__, "..", "data", "life history.bib"), check = :none)

    key_to_citation = map(keys(lifehistorybib), values(lifehistorybib)) do key, bibentry
        authors = bibentry.authors
        auth = if isempty(authors)
            "NO AUTHORS!"
        elseif length(authors) == 1
            authors[1].last
        elseif length(authors) == 2
            authors[1].last * " & " * authors[2].last
        else
            authors[1].last * " et al."
        end
        auth = replace(auth, '{'  => "", '}' => "")

        year = bibentry.date.year
        return key => "$auth ($year)"
    end

    return DataFrame(Key = first.(key_to_citation), Citation = last.(key_to_citation))
end

###  Bayesian models
# A function to couple the data to Bayesian models
# Each of the models themselves are defined just below as functions 
function define_life_history_models(life_history_data)
    # map over hepatica and gigantica data
    map(life_history_data) do data
        # hatching time
        d_ht = data.hatching_time
        # groups ids are by citation and group
        ids = StatsBase.denserank(vcat.(d_ht.Key, d_ht.group))
        model_hatchingtime = censored_gdd_hierarchical(
            d_ht.temperature, d_ht.exp_length, ids, d_ht.observation)

        # prepatent period
        d_pp = data.prepatent_period
        # groups ids are by citation and group
        ids = StatsBase.denserank(vcat.(d_pp.Key, d_pp.group))
        model_prepatentperiod = censored_gdd_hierarchical(
            d_pp.temperature, d_pp.exp_length, ids, d_pp.observation
        )

        # infection efficiency
        # Most sources define infection efficiency as number infected / number exposed
        # one study (for hepatica) uses radioactivity, which is still used
        d_ie = data.infection_efficiency
        d_N = filter(row -> row.unit == "N", d_ie)
        gr_ids = StatsBase.denserank(vcat.(d_N.Key, d_N.group))
        d_rad = filter(row -> startswith("radioactivity")(row.unit), d_ie) # this is empty for gigantica
        model_infectionefficiency = infectionefficiency(
            d_N.temperature, d_N.sample_size, gr_ids,
            d_rad.temperature, d_rad.group,
        ) | (y_bin = d_N.observation, y_rad = d_rad.observation)

        # cercarial release
        d_cr = data.cercarial_release
        ids = StatsBase.denserank(vcat.(d_cr.Key, d_cr.group)) # each citation is a group
        model_cercarialrelease = gaussian_rmax_varying(d_cr.temperature, ids) | 
            (; y = (d_cr.observation))

        # hatching success
        d_hs = data.hatching_success
        gr_ids = StatsBase.denserank(vcat.(d_hs.Key, d_hs.group))
        model_hatchingsuccess = hatching_success(
            d_hs.temperature, d_hs.sample_size, d_hs.observation, gr_ids
        )
        # return all models in a namedtuple
        return (
            hatching_time = model_hatchingtime, 
            prepatent_period = model_prepatentperiod, 
            infection_efficiency = model_infectionefficiency,
            cercarial_release = model_cercarialrelease,
            hatching_success = model_hatchingsuccess
        )
    end
end

#### Actual Bayesian models
# These are the actual models that are used to fit the data
# The @model macro from Turing defines a probabilistic model

# This is used for prepatent period and hatching time
# GDD models but the number of degree days varies per study
import Turing: Normal, LogNormal, MvNormal, Beta, Exponential, Uniform, Binomial, censored,
    arraydist

Turing.@model function censored_gdd_hierarchical(
    T, # temperature of observatoins
    upper, # upper bound of the observations (experiment length)
    ids, # group ids
    y, # observed values
    ngroups = maximum(ids); 
    priors = (
        Tmin = Normal(10, 5),
        days = LogNormal(6, 2)
    )
)
    σ ~ Exponential(10)
    Tmin ~ priors.Tmin
    gdd_mu ~ priors.days
    σ_gdd ~ Exponential(200)
    gdd ~ Turing.filldist(Normal(gdd_mu, σ_gdd), ngroups)
    y_hat = gdd_time.(gdd[ids], Tmin, T; maxval = floatmax())
    for i in eachindex(y)
        if ismissing(upper[i])
            y[i] ~ Normal(y_hat[i], σ)
        else
            y[i] ~ censored(Normal(y_hat[i], σ), nothing, upper[i])
        end
    end
    return FixN(gdd_time, gdd_mu, Tmin)
end

# Infection efficiency model
Turing.@model function infectionefficiency(
    T_bin, n, ids_bin, 
    T_rad, ids_rad, 
    ngroups_rad = maximum(ids_rad; init = 0), ngroups_bin = maximum(ids_bin; init = 0);
    priors = (
        Topt = Normal(20,5), 
        σ = Exponential(10), 
        rmax = Uniform(0.01,0.99)
    )
)
    Topt ~ priors.Topt
    width ~ priors.σ

    rmax_mean ~ priors.rmax
    rmax_sigma ~ censored(Exponential(1), 0, rmax_mean * (1 - rmax_mean) * 0.999)
    x = ((rmax_mean * (1 - rmax_mean) / rmax_sigma) - 1)
    α = rmax_mean * x
    β = (1 - rmax_mean) * x
    rmaxes ~ Turing.filldist(Beta(α, β), ngroups_bin)

    y_hat = gaussian.(rmaxes[ids_bin], Topt, width, T_bin)
    y_bin ~ Turing.arraydist(Binomial.(n, y_hat))

    if ngroups_rad > 0
        rmax_radioactive ~ Turing.filldist(Uniform(0, 1000), ngroups_rad)
        σ ~ Exponential(50)
        y_hat_rad = gaussian.(rmax_radioactive[ids_rad], Topt, width, T_rad)
        y_rad ~ MvNormal(y_hat_rad, σ)
    end

    return FixN(gaussian, rmax_mean, Topt, width)
end

# Used for cercarial release
Turing.@model function gaussian_rmax_varying(
    T, ids, n_group = maximum(ids);
    priors = (
        Topt = Normal(20, 5),
        rmax = LogNormal(6,2),
        breadth = Exponential(5),
        σ_rmax = Exponential(1),
        σ = Exponential(1)
    )
)
    Topt ~ priors.Topt
    rmax_mu ~ priors.rmax
    breadth ~ priors.breadth
    σ_rmax ~ priors.σ_rmax
    σ ~ priors.σ
    rmax ~ Turing.filldist(LogNormal(log(rmax_mu), σ_rmax), n_group)
    y_hat = gaussian.(rmax[ids], Topt, breadth, T)
    y ~ MvNormal(y_hat, σ)
    return FixN(gaussian, rmax_mu, Topt, breadth)
end

# Hatching success
# This one is pretty tricky - it is 0 below some cutoff, then jumps to 1, then 
# Linearly decreases until it starts decreasing exponentially
Turing.@model function hatching_success(
    T, n, k, ids, n_groups = maximum(ids);
    priors = (
        Tmin = Normal(10, 5),
        Tmax = Normal(30, 5),
        k = Exponential(0.5),
        r_Tmax = Uniform(0,1),
        rmax = Uniform(0,1)
    )
)
    Tmin ~ priors.Tmin
    Tmax ~ priors.Tmax
    decl_rate ~ priors.k
    rTmax ~ priors.r_Tmax
    rmax ~ Turing.filldist(priors.rmax, n_groups)
    y_hat = linear_sigmoid.(rTmax, decl_rate, Tmin, Tmax, T; minval = eps()) .* rmax[ids]
    k ~ Turing.product_distribution(Binomial.(n, y_hat))
    return FixN(linear_sigmoid, rTmax, decl_rate, Tmin, Tmax)
end

## Basic functions and utils
gaussian(rmax, Topt, σ, x) = rmax * exp(-(x - Topt)^2 / (2 * σ^2))
function gdd_time(days, Tmin, x; maxval = Inf)
    if x > Tmin
        days / (x - Tmin)
    else
        maxval
    end
end
function linear_sigmoid(r_Tmax, k, Tmin, Tmax, x; minval = zero(x))
    if x < Tmin
        minval
    elseif x > Tmax
        2 * r_Tmax * GLM.logistic(k * (Tmax - x))
    else
        linear(1, (r_Tmax - 1) / (Tmax - Tmin), x - Tmin)
        # 1 - (1 - rate_Tmax) * (x - Tmin) / (Tmax - Tmin)
    end
end
linear(a, b, x) = a + b * x

# FixN struct - to return callable structs with function parameters
struct FixN{F, A}
    f::F
    args::A
end
FixN(f, args...) = FixN(f, args)
(f::FixN)(x) = f.f(f.args..., x)

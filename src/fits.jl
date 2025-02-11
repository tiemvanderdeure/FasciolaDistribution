isvalidunit(unit) = in(unit, ("days", "N", "cercariae/snail")) || Base.contains("radioactivity")(unit)

function load_life_history_data()
    data = CSV.read("data/life_history.csv", DataFrame)
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

### hepatica
function define_life_history_models(life_history_data)
    map(life_history_data) do data
        # hatching time
        d_ht = data.hatching_time
        model_hatchingtime = censored_gdd_hierarchical(
            d_ht.temperature, d_ht.exp_length, StatsBase.denserank(d_ht.Key), d_ht.observation)

        # prepatent period
        d_pp = data.prepatent_period
        model_prepatentperiod = censored_gdd_hierarchical(
            d_pp.temperature, d_pp.exp_length, StatsBase.denserank(d_pp.Key), d_pp.observation
        )

        # infection efficiency
        # a few sources with infection efficiency as %, then one with it expressed in radioactivity units
        d_ie = data.infection_efficiency
        d_N = filter(row -> row.unit == "N", d_ie)
        d_rad = filter(row -> startswith("radioactivity")(row.unit), d_ie)
        groups = StatsBase.denserank(d_rad.unit)
        model_infectionefficiency = infectionefficiency(d_N.temperature, d_N.sample_size, d_N.observation, d_rad.temperature, groups, d_rad.observation)

        # cercarial release
        # tricky to fit a curve, for now just use a gaussian to the log of the cercarial release with varying rmax
        d_cr = data.cercarial_release
        ids = StatsBase.denserank(d_cr.Key)
        model_cercarialrelease = gaussian_rmax_varying(d_cr.temperature, ids) | (; y = (d_cr.observation))

        ## haven't figured this one out yet
        d_hs = data.hatching_success
        model_hatchingsuccess = hatching_success(d_hs.temperature, d_hs.sample_size, d_hs.observation, StatsBase.denserank(d_hs.Key))
        return (
            hatching_time = model_hatchingtime, 
            prepatent_period = model_prepatentperiod, 
            infection_efficiency = model_infectionefficiency,
            cercarial_release = model_cercarialrelease,
            hatching_success = model_hatchingsuccess
        )
    end
        ## or survival
     #   data_survival = filter(row -> row.trait == "survival" && row.units == "N", data_truncatula)
     #   deathrate = -log.(data_survival.observation ./ data_survival.sample_size) ./ data_survival.exp_length


end

function life_history_citations()
    lifehistorybib = parse_file("data/life history.bib", check = :none)
    bibentry = lifehistorybib |> last |> last

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

        year = bibentry.date.year
        return key => "$auth ($year)"
    end

    return DataFrame(Key = first.(key_to_citation), Citation = last.(key_to_citation))
end

#=
function fit_gigantica(data)
    # cercarial release
    # tricky to fit a curve, for now just use a gaussian to the log of the cercarial release with varying rmax
    d_cr = filter(row -> row.trait == "cercarial release" && row.snail == "radix natalensis", data)
    # the radix auricularia data points really mess this one up.
    ids = StatsBase.denserank(d_cr.Key)
    model_cercarialrelease = gaussian_rmax_varying(d_cr.temperature, ids) | (; y = (d_cr.observation))

    ## haven't figured this one out yet - maybe the logistic function is OK?
    d = filter(row -> row.trait == "hatching success", data)
    model_hatchingsuccess = hatching_success(d.temperature, d.sample_size, d.observation, StatsBase.denserank(d.Key))

    ## there isn't really any survival data for radix natalensis...

    return (
        hatchingtime = model_hatchingtime, 
        prepatentperiod = model_prepatentperiod, 
        infectionefficiency = model_infectionefficiency,
        cercarialrelease = model_cercarialrelease,
        hatchingsuccess = model_hatchingsuccess
    )
end
=#

#=
gq = generated_quantities(model_infectionefficiency, chain_infectionefficiency)

model_prepatentperiod = censored_gdd_hierarchical(
    data_prepatent.temperature, data_prepatent.exp_length, StatsBase.denserank(data_prepatent.Key), data_prepatent.observation)
chain_prepatentperiod = sample(model_prepatentperiod, NUTS(), 2000, progress=true);



d = data_hatchingsuccess
model_hathcingsuccess = hatching_success(d.temperature, d.sample_size, d.observation, StatsBase.denserank(d.Key))
chain_hatchingsuccess = sample(model_hathcingsuccess, NUTS(), 1000, progress=true);

plot(exp.(gaussian.(5.8, 24, 18.8, 1:40)))


model_prepatentperiod2 = censored_gdd(
    data_prepatent.temperature, data_prepatent.exp_length, data_prepatent.observation;
    gdd_prior = LogNormal(6,2))
chain_prepatentperiod2 = sample(model_prepatentperiod2, NUTS(), 1000, progress=true);

model_hatchingtime2 = censored_gdd(
    data_hatchingtime.temperature, data_hatchingtime.exp_length, data_hatchingtime.observation)
chain_hatchingtime2 = sample(model_hatchingtime2, NUTS(), 1000, progress=true);

model_survival2 = survival(d.temperature, d.exp_length, d.sample_size, d.observation)
chain_survival = sample(model_survival2, NUTS(), 1000, progress=true);

#=
for (data, fasciola) in zip([data_gig_radix, data_hep_galba], ["gigantica", "hepatica"])
    for trait in ["prepatent period", "infection efficiency", "cercarial release", "hatching success", "hatching time",]
        traitdata = filter(row -> row.trait == trait, data)
        if nrow(traitdata) > 0
            data_transformed = map(eachrow(traitdata)) do row
                row.units == row.units == "N" ? row.observation ./ row.sample_size : row.observation
            end 
            citationkeys = unique(traitdata.Key)
            nkeys = length(citationkeys)
            ylabel = traitdata.units[1] == "N" ? "share" : traitdata.units[1]

            fig = Figure()
            ax = Axis(
                fig[1,1]; 
                title = trait, xlabel = "temperature (°C)", ylabel, 
                limits = (nothing,nothing,0,nothing))
            for (i, key) in enumerate(citationkeys) 
                indices = findall(traitdata.Key .== key)
                plotdata = traitdata[indices, :]
                marker = if trait in ["prepatent period", "hatching time"]
                    ifelse.(.~ismissing.(plotdata.exp_length) .&& plotdata.exp_length .<= plotdata.observation, :utriangle, :circle) 
                else  
                    :circle
                end
                scatter!(
                    ax, plotdata.temperature, data_transformed[indices]; 
                    color = i, colorrange = (1, nkeys), colormap = :rainbow,
                    marker, label = key)
            end
            fig[2,1] = Legend(fig, ax; tellheight = true, tellwidth = false)

            save("plots/$fasciola/$trait.png", fig)
        end
    end
end
=#
=#


#=
scatter(data_cercarialrelease.temperature, log.(data_cercarialrelease.observation); color = collect(ids))

scatter(data_infectionefficiency.temperature, data_infectionefficiency.observation ./ data_infectionefficiency[!, "sample size"] )

lambda = log(2)/80
1/lambda
exp(-600*lambda)



model_prepatent = censored_gdd(data_prepatent.temperature, data_prepatent.exp_length, data_prepatent.observation)
chain_prepatent = sample(model_prepatent, sampler, 1000, progress=true);

model_infectionefficiency = linear_model_binomial(data_infectionefficiency.temperature, data_infectionefficiency[!, "sample size"], data_infectionefficiency.observation)
chain_infectionefficiency = sample(model_infectionefficiency, sampler, 1000, progress=true);


T = 1.0:40.0
chain = chain_infectionefficiency
pred = GLM.logistic.(chain[:intercept] .+ chain[:a] .* T')
pred_mean = vec(median(pred, dims=1))

fig, ax, sc = scatter(data_infectionefficiency.temperature, data_infectionefficiency.observation ./  data_infectionefficiency[!, "sample size"])
lines!(ax, pred_mean)
fig


data_prepatent = copy(filter(row -> row.trait == "prepatent period", data_gig))
data_prepatent_hep = copy(filter(row -> row.trait == "prepatent period", data_hep))
data_hatchingtime_hep = copy(filter(row -> row.trait == "hatching time", data_hep))
data_surv_natalensis = copy(filter(row -> row.trait == "survival", data_incl))

data_cercarialrelease = copy(filter(row -> row.trait == "cercarial release", data_gig))
scatter(data_cercarialrelease.temperature, data_cercarialrelease.observation; color = Vector(StatsBase.denserank(data_cercarialrelease.Key)))

sampler = NUTS()

model_prepatent_gig = censored_gdd(data_prepatent.temperature, data_prepatent.exp_length, data_prepatent.observation)
chain_prepatent_gig = sample(model_prepatent_gig, sampler, 1000, progress=true);

model_prepatent_hep = censored_gdd(data_prepatent_hep.temperature, data_prepatent_hep.exp_length, data_prepatent_hep.observation)
chain_prepatent_hep = sample(model_prepatent_hep, sampler, 1000, progress=true);


# calculate the expected y for each temperature between 10 and 40 C 
T = 1.0:40.0
pred = gdd_time.(vec(chain[:Tmin]), vec(chain[:gdd]), T')
pred_mean = vec(median(pred, dims=1))

fig, ax, sc = scatter(data_prepatent.temperature, data_prepatent.observation);
limits!(ax, 0, 35, 0, 250)
lines!(ax, T[12:end], pred_mean[12:end])
fig


pred = gdd_time.(vec(chain_prepatent_hep[:Tmin]), vec(chain_prepatent_hep[:gdd]), T')
pred_mean = vec(median(pred, dims=1))
marker = ifelse.(.~ismissing.(data_prepatent_hep.exp_length) .&& data_prepatent_hep.exp_length .<= data_prepatent_hep.observation, :utriangle, :circle)
color = Vector(StatsBase.ordinalrank(data_prepatent_hep.Key))
fig, ax, sc = scatter(data_prepatent_hep.temperature, data_prepatent_hep.observation; color, marker)
limits!(ax, 5, 42, 0, 250)
i = findfirst(pred_mean .< 1000)
lines!(ax, T[i:end], pred_mean[i:end])
fig


scatter(data_hatchingtime_hep.temperature, data_hatchingtime_hep.observation)
=#
using GLMakie, FasciolaDistribution, Rasters

global figurepath = joinpath("images")
global supplementalspath = joinpath("images", "supplementary")

const FD = FasciolaDistribution # any function called with FD. is from this package

global snails = (:galba, :radix)
global fasc_sps = (:hepatica, :gigantica)
global temperatures = 5.0:0.1:40.0 # temperatures to plot


#### Figure 1 - occurrence data
fig1 = let
    samplingbgcoordinates = values.(mapreduce(x -> getfield(x, :geo), vcat, sbg)),
    radixcoordiantes = values.(occs.radix.geo),
    galbacoordinates = [values.(occs.galba_eu.geo); values.(occs.galba_af.geo)]

    fig = Figure(size = (1000, 700))
    ax = Axis(
        fig[1, 1], 
        limits = (-25, 53, -35, 73), 
        aspect = AxisAspect(1), 
        title = "Recorded occurrences of Lymnaeids"
    ) # probably use geoaxis instead
    hidedecorations!(ax); hidespines!(ax)
    poly!(ax, countries.geometry, color = :transparent, strokewidth = 0.7)

    scatter!(samplingbgcoordinates; color = :grey, label = "Other Lymnaeids", markersize = 6)
    scatter!(radixcoordiantes; color = :red, label = rich("Radix natalensis", font = :italic), markersize = 6)
    scatter!(galbacoordinates; color = :blue, label = rich("Galba truncatula", font = :italic), markersize = 6)
    axislegend(ax, position = (0.15, 0.15))
    save(joinpath(figurepath, "figure1.png"), fig)
    fig
end

### Figure 2 - life history posteriors
as_label(x::Symbol) = replace(string(x), "_" => " ")

# plot for these temperatures
fig2 = let temperatures = 5.0:0.1:40.0,
    colors = (hepatica = :blue, gigantica = :red),
    limits = (hatching_time = (0, 100), infection_efficiency = (0, 1), prepatent_period = (0, 200), cercarial_release = (0, 1500), hatching_success = (0, 1)),
    ylabels = (hatching_time = "days", infection_efficiency = "share", prepatent_period = "days", cercarial_release = "cercariae", hatching_success = "share")

    fig = Figure(fontsize = 12, size = (850, 450))

    ncol = 3

    for (i, trait) in enumerate(trait_keys)
        ax = Axis(fig[(i-1) ÷ ncol + 1, mod1(i, ncol)]; 
            limits = (extrema(temperatures)..., limits[trait]...),
            title = as_label(trait),
            xlabel = "temperature (°C)", ylabel = ylabels[trait])
        for fasciola in keys(life_history_data)
            preds = temperatures' .|> gqs[fasciola][trait]
            plot_quantiles!(ax, temperatures, preds, color = colors[fasciola])
        end
    end
    for fasciola in keys(life_history_data)
        preds = life_cycle_model.(temperatures', gqs_rt[fasciola])
        # normalize these to be on a 0-1 scale
        preds_normalized = preds ./ maximum(preds, dims = 2)
        # plot and save
        ax = Axis(fig[2, 3], limits = (extrema(temperatures)..., 0, 1), 
            title = "transmission strength", xlabel = "temperature (°C)", ylabel = "cercariae per egg (rel.)")
        ls = plot_quantiles!(ax, temperatures, preds_normalized, color = colors[fasciola])
    end
    for i in 1:6 
        Label(fig[(i-1) ÷ ncol + 1, mod1(i, ncol), TopLeft()], "($(Char(i+96)))", font = :bold,
            tellheight = false, tellwidth = false, halign = :left, valign = :top)
    end

    # Legend
    le_f = [LineElement(; color) for color in colors]
    f_labels = [rich("F. $f"; font = :italic) for f in keys(colors)]
    le_post = [LineElement(; linestyle) for linestyle in [:solid, :dash]]
    post_labels = ["mean", "95% CI"]
    Legend(fig[1:2,ncol+1], [le_f, le_post], [f_labels, post_labels], ["Species", "Posterior"])

    fig
end
save(joinpath(figurepath, "figure_2.png"), fig2)


### Figure 3 - current 
hydro_all, host_all, trans_all, temp_all, risk_all =
    FasciolaDistribution.read_predictions((gcms, ghms, dates, ssps); lazy = true)

hydro, host, trans, temp, risk =
    FasciolaDistribution.read_predictions(; prefix = "mean_")
ext =  extent(risk.current.hepatica)
ext = Extent((X = (ext.X[1] - 0.00001, ext.X[2] + 0.00001)), Y = (ext.Y[1] - 0.00001, ext.Y[2] + 0.00001))
asp = (ext.X[2]- ext.X[1]) / (ext.Y[2] - ext.Y[1]) # aspect ratio of maps!
    
_risk = map(temp, host) do t, h
    map(values(t), values(h)) do t, h
        sqrt.(t .* h)
    end |> NamedTuple{keys(t)}
end

# current conditions

fig3 = let fig = Figure(fontsize = 14, size = (1100, 600))
    plot_kw = (
        colorrange = (0.001,1),
        colormap = Reverse(:Spectral),
        lowclip = :black
    )

    for i in 1:2
        snail = snails[i]
        fasc_sp = fasc_sps[i]

        ax1 = Axis(fig[i, 1])
        plot!(ax1, hydro.current; plot_kw...)
        ax2 = Axis(fig[i, 2])
        plot!(ax2, temp.current[fasc_sp]; plot_kw...)
        ax3 = Axis(fig[i, 3])
        plot!(ax3, trans.current[fasc_sp]; plot_kw...)
        ax4 = Axis(fig[i, 4])
        plot!(ax4, host.current[snail]; plot_kw...)
        ax5 = Axis(fig[i, 5])
        plot!(ax5, risk.current[fasc_sp]; plot_kw...)

        for ax in (ax1, ax2, ax3, ax4, ax5)
            hidedecorations!(ax); hidespines!(ax)
        end
        
        Label(fig[i, 0], "Fasciola " * string(fasc_sp); tellheight = false,rotation = pi/2, font = :bold_italic)
        rowsize!(fig.layout, 1, Aspect(1, 1/asp))
    end

    labels = ["Hydrological suitability", "Temperature suitability", "Parasite suitability", "Host suitability", "Transmission risk"]
    for i in 1:5
        Label(fig[0, i], labels[i]; tellwidth = false, font = :bold)
    end

    Colorbar(fig[1:2, 6]; plot_kw...)
    resize_to_layout!(fig)

    save(joinpath(figurepath, "figure3.png"), fig)
    fig
end

#### Alternative figure 3
fig3 = let fig = Figure(fontsize = 14, size = (600, 600))
    snails = (:galba, :radix)
    fasc_sps = (:hepatica, :gigantica)

    plot_kw = (
        colorrange = (0,1),
        colormap = Reverse(:Spectral),
    )

    for i in 1:2
        snail = snails[i]
        fasc_sp = fasc_sps[i]

        ax1 = Axis(fig[i, 1])
        plot!(ax1, temp.current[fasc_sp]; plot_kw...)
        ax2 = Axis(fig[i, 2])
        plot!(ax2, host.current[snail]; plot_kw...)
        ax3 = Axis(fig[i, 3])
        plot!(ax3, _risk.current[fasc_sp]; plot_kw...)

        for ax in (ax1, ax2, ax3)
            hidedecorations!(ax); hidespines!(ax)
        end
        
        Label(fig[i, 0], "Fasciola " * string(fasc_sp); tellheight = false,rotation = pi/2, font = :bold_italic)
        rowsize!(fig.layout, i, Aspect(1, 1/asp))
    end

    labels = ["Temperature suitability", "Host suitability", "Transmission risk"]
    for i in 1:3
        Label(fig[0, i], labels[i]; tellwidth = false, font = :bold)
    end

    Colorbar(fig[1:2, 4]; plot_kw...)
    resize_to_layout!(fig)

    save(joinpath(figurepath, "figure3_alt.png"), fig)
    fig
end

### Future predictions
# generates figures for temperature, host suitability and risk
# under SSP126 and SSP370
# risk for SSP370 is in the main manuscript, all others are in the supplementals
(figs3, figs4), (figs5, figs6), (figs2, fig4) = let 
    data = (temp, host, _risk)
    legend_labels = ("Temperature suitability", "Host suitability", "Transmission risk")
    fasc_labels = ("F. hepatica", "F. gigantica")
    snail_labels = ("Galba truncatula", "Radix natalensis")
    x_labels = (fasc_labels, snail_labels, fasc_labels)
    map(data, legend_labels, x_labels) do data, label, species
        map(ssps) do ssp
            fig = Figure(fontsize = 14, size = (800, 600))
            kws = (
                colorrange = (0,1),
                colormap = Reverse(:Spectral)
            )

            for s in 1:2
                fasciola_sp = fasc_sps[s]

                ax_current = Axis(fig[s, 1])
                plot!(ax_current, data.current[s]; kws...)
                hidedecorations!(ax_current); hidespines!(ax_current)

                for d in eachindex(dates)
                    ax = Axis(fig[s, d+1])
                    plot!(ax, data.future[s][ssp = At(ssp), Ti = d]; kws...)
                    hidedecorations!(ax); hidespines!(ax)
                    if s == 1
                        d1, d2 = val(Rasters.span(dates))[d, :]
                        Label(fig[0, d+1], "$(year(d1))-$(year(d2))"; tellwidth = false, font = :bold)
                    end
                end
                Label(fig[s, 0], species[s]; tellheight = false,rotation = pi/2, font = :bold_italic)
            end
            Label(fig[0, 1], "current"; tellwidth = false, font = :bold)

            Colorbar(fig[1:2, 4]; kws..., label)

            fig
        end
    end
end
save(joinpath(figurepath, "figure4_alt.png"), fig4)

save(joinpath(supplementalspath, "forecast_ssp126.png"), figs2)
save(joinpath(supplementalspath, "temp_ssp126.png"), figs3)
save(joinpath(supplementalspath, "temp_ssp370.png"), figs4)
save(joinpath(supplementalspath, "host_ssp126.png"), figs5)
save(joinpath(supplementalspath, "host_ssp370.png"), figs6)


### Overlap of fasciola risk and livestock
sheep = Raster("D:/data/FAO/GLW4-2020.D-DA.SHP.tif")
cattle = Raster("D:/data/FAO/GLW4-2020.D-DA.CTL.tif")
sheep = disaggregate(crop(sheep; to = ext), 2) # to match resolution of risk
cattle = disaggregate(crop(cattle; to = ext), 2) # to match resolution of risk

livestock = (; sheep, cattle)

fig5 = let fig = Figure(fontsize = 10)

    xticks = [0, 0.25, 0.5, 0.75, 1]
    yticks = [0, 10, 20, 50, Inf]
    tolochko_redblue = [
        colorant"#dd0027" colorant"#4f2d4c"
        colorant"#dcdcdc" colorant"#0072a9"
    ]
    cmap = FD.bivariate_colormap(xticks, yticks; colors = tolochko_redblue)

    for ti in 1:2
        ti_idx = (ti-1)*2
        for f_idx in 1:2
            fs = fasc_sps[f_idx]
            for l_idx = 1:2
                ax = Axis(fig[l_idx, f_idx+ti_idx])
                l_sp = keys(livestock)[l_idx]
                plot!(
                    ax,
                    FD.rasters_to_cmap(_risk[ti][fs][Ti = 2, ssp = 2], livestock[l_sp]; cmap)
                )           
                hidedecorations!(ax); hidespines!(ax) 
                if ti == 1 && f_idx == 1
                    Label(fig[l_idx, 1, Left()], string(l_sp); tellheight = false, rotation = pi/2, font = :bold)
                end
            end
            Label(fig[0, f_idx+ti_idx], "F. " * string(fs); tellwidth = false, font = :bold_italic)
        end
        Box(fig[-1, (ti*2-1):ti*2], color = (:lightgrey, 0.5), strokevisible = false)
        timelabel = ti == 1 ? "Current" : "Future"
        Label(fig[-1, (ti*2-1):ti*2], timelabel; tellwidth = false, font = :bold)
    end
    for i in 1:2
            # aspect ratio
            rowsize!(fig.layout, i, Aspect(1, 1/asp))
    end
    FD.bivariate_cmap_legend(
        fig[1:2, 5], cmap; 
        aspect = AxisAspect(1),
        yticksf = x -> isfinite(x) ? string(Int(x)) : ">$(Int(maximum(filter(isfinite, yticks))))",
        ylabel = "Livestock density", xlabel = rich(rich("Fasciola", font = :italic), " risk"),
        xticklabelrotation = pi/4,
        alignmode = Mixed(bottom = 0) # hack to include protrusions so we don't get excess white space
    )

    colsize!(fig.layout, 5, Relative(0.15))

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    resize_to_layout!(fig)
    fig 
end
resize_to_layout!(fig5)
save(joinpath(figurepath, "risk_livestock.png"), fig5)


#### Supplementals for life history traits
### plot life history trait data and posterior estimates of curves against temperature
using DataFrames: leftjoin

citations = FasciolaDistribution.life_history_citations()

# Sample from the models. Priors are identical so only need to sample from 1
life_history_priors = map(life_history_models.hepatica) do m
    chain = sample(Xoshiro(0), m, Prior(), 10_000, progress=true)
end
gqs_prior = map(life_history_models.hepatica, life_history_priors) do m, c
    vec(generated_quantities(m, c))
end

for fasciola in fasc_sps
    dir = joinpath(supplementalspath, String(fasciola))
    isdir(dir) || mkdir(dir)

    # loop through traits
    #(trait, trait_str) = first(zip(trait_keys, traits_str))
    for (trait, trait_str) in zip(trait_keys, traits_str)
        @show trait
        fig = Figure()
        ax = Axis(
            fig[1,1]; 
            title = rich(rich("F. $(fasciola)"; font = :bold_italic),  " $(trait_str)"),
            alignmode = Mixed(bottom = 0))
        legend_gl = GridLayout(fig[1,2])

        # get and plot prior distributions
        preds_prior = temperatures' .|> gqs_prior[trait]
        pr1 = FD.shade_quantiles!(ax, temperatures, preds_prior; p = [0.25,0.75], color = (:gray, 0.4))
        pr2 = FD.shade_quantiles!(ax, temperatures, preds_prior; p = [0.025,0.975], color = (:gray, 0.2))
        
        # get and plot the posterior predictions
        preds = temperatures' .|> gqs[fasciola][trait]
        ls = plot_quantiles!(ax, temperatures, preds)
        # get and plot the data
        data = life_history_data[fasciola][trait]
        data_w_cit = leftjoin(data, citations, on = :Key)
        plot_life_history!(
            ax, data_w_cit, temperatures; 
            legend_gp = legend_gl[1,1]
        )
        Legend(legend_gl[2,1], [pr1, pr2], ["50% CI", "95% CI"]; tellheight = true)
        Legend(legend_gl[3,1], ls[[3,2]], ["mean", "95% CI"]; tellheight = true)
        Label(legend_gl[1,1, Top()], "Data sources"; font = :bold)
        Label(legend_gl[2,1, Top()], "Prior"; font = :bold)
        Label(legend_gl[3,1, Top()], "Posterior"; font = :bold)

        colsize!(fig.layout, 1, Fixed(450))
        rowsize!(fig.layout, 1, Fixed(450))
        resize_to_layout!(fig)
        # save with posterior
        save(joinpath(supplementalspath, string(fasciola), "$trait.png"), fig)
    end
end

## sensitivity analysis to assumptions on egg and snail death rate
fasciola = :hepatica
snail_death_rates = (0.002, 0.01, 0.05)
egg_death_rates = (0.002, 0.01, 0.05)
colors = (hepatica = :blue, gigantica = :red)

fig_sensitivity_analysis = let fig = Figure(size = (800, 600))
    for (i, snail_death_rate) in enumerate(snail_death_rates)
        for (j, egg_death_rate) in enumerate(egg_death_rates)
            isbottom = i == 3
            isleft = j == 1
            xlabel = isbottom ?  "temperature (°C)" : ""
            ylabel = isleft ? "cercariae per egg (rel.)" : ""
            ax = Axis(
                fig[i, j]; 
                limits = (extrema(temperatures)..., 0, 1), 
                xlabel, ylabel,
                xticklabelsvisible = isbottom, yticklabelsvisible = isleft
            )

            for fasciola in fasc_sps
                preds = life_cycle_model.(
                    temperatures', gqs_rt[fasciola]; 
                    snail_death_rate, egg_death_rate
                )
                # normalize these to be on a 0-1 scale
                preds_normalized = preds ./ maximum(preds, dims = 2)
                # plot and save
                ls = plot_quantiles!(ax, temperatures, preds_normalized, color = colors[fasciola])
            end
            if i == 1
                Label(fig[0, j], "egg death rate: $(egg_death_rate)"; tellwidth = false, font = :bold)
            end
        end
        Label(fig[i, 0], "snail death rate: $(snail_death_rate)"; tellheight = false,rotation = pi/2, font = :bold)
    end
    fig
end
save(joinpath(supplementalspath, "sensitivity_analysis.png"), fig_sensitivity_analysis)


#=
fig = Figure()
Label(
    fig[0, 1:4], 
    rich("Climatic suitability for ", rich("Fasciola", font = :bold_italic), " transmission"); 
    font = :bold, fontsize = 16)
season_labels = ["Dec-Feb", "Mar-May", "Jun-Aug", "Sep-Nov"]
colorrange = (0,1); cmap = Reverse(:Spectral)
maps = GridLayout(fig[1:2, 1:4])

for (i, fasciola) in enumerate(keys(transmission_quarterly))
    transmission_raster_quarterly = transmission_quarterly[fasciola]
    Label(fig[i, 0], string(fasciola); tellheight = false,rotation = pi/2, font = :bold_italic)
    for q in 1:4
        ax = Axis(maps[i, q])
        plot!(ax, @view(transmission_raster_quarterly[Dim{:quarter}(q)]); colorrange, colormap = cmap)
        plot!(ax, desert_mask; colormap = [:lightgrey, :transparent])
        hidedecorations!(ax); hidespines!(ax)
        if i == 1
            Label(fig[1, q, Top()], season_labels[q]; tellwidth = false, font = :bold)
        end
    end
end
rowgap!(maps, 0); colgap!(maps, 0)
fig
save("images/quarterly_transmission.png", fig)


# for on the website
fasciola = :hepatica
host = :galba
times = [:current, :future]
countries_af = get_countries()

for (fasciola, host) in zip(keys(life_history_data), keys(occurrences))

    fig = Figure()

    maps = GridLayout(fig[1, 1])

    colorrange = (0,1); cmap = Reverse(:Spectral)
    ext = Rasters.Extents.extent(bio.current)
    limits = (ext.X..., ext.Y...)

    for (i, time) in enumerate(times)
        Label(maps[i, 0], string(time); tellheight = false,rotation = pi/2, font = :bold)
        ax = Axis(maps[i,1]; limits, title = ifelse(i == 1, "host suitability", ""))
        plot!(ax, host_suitability[time][host]; colorrange, colormap = cmap)

        ax2 = Axis(maps[i,2]; limits, title = ifelse(i == 1, "transmission strength", ""))
        plot!(ax2, temperature_suitability_annual[time][fasciola]; colorrange, colormap = cmap)

        ax3 = Axis(maps[i,3]; limits, title = ifelse(i == 1, "infection risk", ""))
        plot!(ax3, temperature_suitability_annual[time][fasciola] .* host_suitability[time][host]; colorrange, colormap = cmap)

        for ax in (ax, ax2, ax3)
            hidespines!(ax); hidedecorations!(ax)
            poly!(ax, countries_af, color = :transparent, strokewidth = 0.3)
        end
    end

    Colorbar(maps[1:2,4]; colorrange, colormap = cmap)

    save("images/infectionrisk_website_$fasciola.png", fig)
end

## Current climate, quarterly transmission
fig = Figure()
Label(
    fig[0, 1:4], 
    rich("Climatic suitability for ", rich("Fasciola", font = :bold_italic), " transmission"); 
    font = :bold, fontsize = 16)
season_labels = ["Dec-Feb", "Mar-May", "Jun-Aug", "Sep-Nov"]
colorrange = (0,1); cmap = Reverse(:Spectral)
maps = GridLayout(fig[1:2, 1:4])

for (i, fasciola) in enumerate(keys(transmission_quarterly))
    transmission_raster_quarterly = transmission_quarterly[fasciola]
    Label(fig[i, 0], string(fasciola); tellheight = false,rotation = pi/2, font = :bold_italic)
    for q in 1:4
        ax = Axis(maps[i, q])
        plot!(ax, @view(transmission_raster_quarterly[Dim{:quarter}(q)]); colorrange, colormap = cmap)
        plot!(ax, desert_mask; colormap = [:lightgrey, :transparent])
        hidedecorations!(ax); hidespines!(ax)
        if i == 1
            Label(fig[1, q, Top()], season_labels[q]; tellwidth = false, font = :bold)
        end
    end
end
rowgap!(maps, 0); colgap!(maps, 0)
fig
save("images/quarterly_transmission.png", fig)


## Future climate, hepatica, quarterly transmission
fig = Figure()
Label(
    fig[0, 1:4], 
    rich("Climatic suitability for ", rich("Fasciola hepatica", font = :bold_italic), " transmission"); 
    font = :bold, fontsize = 16)
season_labels = ["Dec-Feb", "Mar-May", "Jun-Aug", "Sep-Nov"]
colorrange = (0,1); cmap = Reverse(:Spectral)
maps = GridLayout(fig[1:2, 1:4])

time = ["current", "2081-2100"]
for (i, transmission_raster_quarterly) in enumerate((transmission_quarterly.hepatica, transmission_quarterly_future.hepatica))
    Label(fig[i, 0], time[i]; tellheight = false,rotation = pi/2, font = :bold)
    for q in 1:4
        ax = Axis(maps[i, q])
        plot!(ax, @view(transmission_raster_quarterly[Dim{:quarter}(q)]); colorrange, colormap = cmap)
        plot!(ax, desert_mask; colormap = [:lightgrey, :transparent])
        hidedecorations!(ax); hidespines!(ax)
        if i == 1
            Label(fig[1, q, Top()], season_labels[q]; tellwidth = false, font = :bold)
        end
    end
end
rowgap!(maps, 10); colgap!(maps, 10)
fig
save("images/quarterly_transmission_hepatica_future.png", fig)

=#


@model function loop_bern()
    y = Vector{Bool}(undef, 10)
    p ~ Uniform(0, 1)
    for i in 1:10
        y[i] ~ Bernoulli(p)
    end
end

model = loop_bern() | (; y = trues(10))
sample(model, NUTS(), 1000)

@model function demo_mv(::Type{TV}=Float64) where {TV}
    m = Vector{TV}(undef, 2)
    m .~ Normal.([1.0,1.0], [1.0,1.0])
    return m
end

model = demo_mv();

conditioned_model = condition(model, m = [2.0, 1.0]);
conditioned_model = model | (; m = [missing, 2.0])
sample(conditioned_model, NUTS(), 10)

@model function demo_mv(m)
 #   m = Vector{Float64}(undef, 2)
    m[1] ~ Normal(0, 1)
    m[2] ~ Normal(0, 1)
end

model = demo_mv(rand(2))# | (; m = [1, 2])
sample(model, NUTS(), 1000)
sample(demo_mv() | m = [1, 2], NUTS(), 1000)
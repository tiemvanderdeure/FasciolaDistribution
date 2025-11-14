using CairoMakie, FasciolaDistribution, Rasters, RasterDataSources, Statistics, NaturalEarth
import Dates: year
import GeometryBasics
const FD = FasciolaDistribution
ENV["FASCIOLA_DATA_PATH"] = joinpath(RasterDataSources.rasterpath(), "FasciolaDistribution")

global figurepath = joinpath("images")
global supplementalspath = joinpath("images", "supplementary")

global snails = (:galba, :radix)
global fasc_sps = (:hepatica, :gigantica)
global temperatures = 5.0:0.1:40.0 # temperatures to plot

hydro_all, host_all, trans_all, temp_all, risk_all =
    FD.read_predictions((FD.gcms, FD.ghms, FD.dates, FD.ssps); lazy = true)

# load mean predictions, and get rid of a silly corner of Greenland
greenland = GeometryBasics.Polygon(GeometryBasics.Rect2f((-25, 67), (10, 10)))

meanpreds = FD.read_predictions(; prefix = "mean_")

hydro = map(x -> mask(x; with = greenland, invert = true), meanpreds[1])
host, trans, temp, risk = 
    map(x -> mapmap(x2 -> mask(x2; with = greenland, invert = true), x), meanpreds[2:end])

FD.mapmap(risk) do x
    x[.!ismissing.(x) .&& isnan.(x)] .= missing
end

ext = extent(risk.current.hepatica)
global asp = (ext.X[2]- ext.X[1]) / (ext.Y[2] - ext.Y[1]) # aspect ratio of maps!
global px_per_unit = 5 # for high-res images

### Figure 1 - life history posteriors
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
            title = "transmission capacity", xlabel = "temperature (°C)", ylabel = "cercariae per egg (rel.)")
        ls = plot_quantiles!(ax, temperatures, preds_normalized, color = colors[fasciola])
    end
    for i in 1:6 
        FD.Label_subplot(fig[(i-1) ÷ ncol + 1, mod1(i, ncol), TopLeft()]; n = i, valign = :top)
    end

    # Legend
    le_f = [LineElement(; color) for color in colors]
    f_labels = [rich("F. $f"; font = :italic) for f in keys(colors)]
    le_post = [LineElement(; linestyle) for linestyle in [:solid, :dash]]
    post_labels = ["mean", "95% CI"]
    Legend(fig[1:2,ncol+1], [le_f, le_post], [f_labels, post_labels], ["Species", "Posterior"])

    fig
end
save(joinpath(figurepath, "figure2.png"), fig2; px_per_unit)

### Figure 2 - current 
# current conditions

fig3 = let 
    fig = Figure(fontsize = 14, size = (930, 770))
    plot_kw = (
        colorrange = (0,1),
        colormap = Reverse(:Spectral),
    )

    components = GridLayout(fig[1,1])
    riskindex = GridLayout(fig[2,1])

    titlegap = 7

    # plot all the components
    ax_temp_hep = Axis(components[1, 1]; titlegap,
        title = rich(rich("F. hepatica", font = :bold_italic), "\ntransmission capacity"), titlelineheight = 0)
    plot!(ax_temp_hep, temp.current.hepatica; plot_kw...)
    ax_host_hep = Axis(components[1, 2]; titlegap,
        title = rich(rich("G. truncatula", font = :bold_italic), "\nclimatic suitability"))
    plot!(ax_host_hep, host.current.galba; plot_kw...)
    ax_hydro = Axis(components[1, 3]; titlegap, title = "Hydrological suitability")
    plot!(ax_hydro, hydro.current; plot_kw...)
    ax_temp_gig = Axis(components[1, 4]; titlegap,
        title = rich(rich("F. gigantica", font = :bold_italic), "\ntransmission capacity"))
    plot!(ax_temp_gig, temp.current.gigantica; plot_kw...)
    ax_host_gig = Axis(components[1, 5]; titlegap,
        title = rich(rich("R. natalensis", font = :bold_italic), "\nclimatic suitability"))
    plot!(ax_host_gig, host.current.radix; plot_kw...)

    # plot the actual risk indices
    ax_risk_hep= Axis(riskindex[1, 1:3]; aspect= asp, titlesize = 16,
        title = rich(rich("F. hepatica", font = :bold_italic), "\ntransmission risk"))
    plot!(ax_risk_hep, risk.current.hepatica; plot_kw...)
    ax_risk_gig = Axis(riskindex[1, 3:5]; aspect= asp, titlesize = 16,
        title = rich(rich("F. gigantica", font = :bold_italic), "\ntransmission risk"))
    plot!(ax_risk_gig, risk.current.gigantica; plot_kw...)
    
    for ax in fig.content
        hidedecorations!(ax); hidespines!(ax)
    end
    
    for i in 1:5 FD.Label_subplot(components[1,i]; n= i, valign = 1.04) end
    FD.Label_subplot(riskindex[1,1:3]; n=6, halign = 1/4)
    FD.Label_subplot(riskindex[1,3:5]; n=7, halign = 1/4)
    rowsize!(components, 1, Aspect(1, 1/asp))
    
    Colorbar(fig[1:2, 2]; height = Relative(0.7), valign = 0.4, plot_kw...)

    Label(fig[1, 0], "Risk index components", rotation = pi/2, 
        tellheight = false, font = :bold, fontsize = 16)
    Label(fig[2, 0], "Integrated risk index", rotation = pi/2, 
        tellheight = false, font = :bold, fontsize = 16)

    resize_to_layout!(fig)

    fig
end
save(joinpath(figurepath, "figure3_noarrows.png"), fig3; px_per_unit) # add arrows manually

### Future predictions
# generates figures for temperature, host suitability and risk
# under SSP126 and SSP370
# risk for SSP370 is in the main manuscript, all others are in the supplementals
(figs3, figs4), (figs5, figs6), (figs2, fig4) = let 
    data = (temp, host, risk)
    legend_labels = ("Temperature suitability", "Snail host suitability", "Transmission risk")
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
                FD.Label_subplot(fig[s,1], ((s-1)*3+1))

                for d in eachindex(dates)
                    ax = Axis(fig[s, d+1])
                    plot!(ax, data.future[s][ssp = At(ssp), Ti = d]; kws...)
                    hidedecorations!(ax); hidespines!(ax)
                    if s == 1
                        d1, d2 = val(Rasters.span(dates))[d, :]
                        Label(fig[0, d+1], "$(year(d1))-$(year(d2))"; tellwidth = false, font = :bold)
                    end
                    FD.Label_subplot(fig[s,d+1]; n = (s-1)*3+d+1)
                end
                Label(fig[s, 0], species[s]; tellheight = false,rotation = pi/2, font = :bold_italic)
                rowsize!(fig.layout, s, Aspect(1, 1/asp))
                resize_to_layout!(fig)
            end
            Label(fig[0, 1], "current"; tellwidth = false, font = :bold)

            Colorbar(fig[1:2, 4]; kws..., label)

            fig
        end
    end
end;
save(joinpath(figurepath, "figure4.png"), fig4; px_per_unit)

save(joinpath(supplementalspath, "forecast_ssp126.png"), figs2; px_per_unit)
save(joinpath(supplementalspath, "temp_ssp126.png"), figs3; px_per_unit)
save(joinpath(supplementalspath, "temp_ssp370.png"), figs4; px_per_unit)
save(joinpath(supplementalspath, "host_ssp126.png"), figs5; px_per_unit)
save(joinpath(supplementalspath, "host_ssp370.png"), figs6; px_per_unit)

## Supplemental with hydrological suitability
for ssp in ssps
    fig = Figure(fontsize = 14, size = (800, 600))
    kws = (
        colorrange = (0,1),
        colormap = Reverse(:Spectral)
    )

    ax_current = Axis(fig[1, 1])
    plot!(ax_current, hydro.current; kws...)
    hidedecorations!(ax_current); hidespines!(ax_current)

    for d in eachindex(dates)
        ax = Axis(fig[1, d+1])
        plot!(ax, hydro.future[ssp = At(ssp), Ti = d] .- hydro.current; kws...)
        hidedecorations!(ax); hidespines!(ax)
        d1, d2 = val(Rasters.span(dates))[d, :]
        Label(fig[0, d+1], "$(year(d1))-$(year(d2))"; tellwidth = false, font = :bold)
    end
    Label(fig[0, 1], "current"; tellwidth = false, font = :bold)
    Colorbar(fig[1, 4]; kws..., label = "Hydrological suitability")

    rowsize!(fig.layout, 1, Aspect(1, 1/asp))
    resize_to_layout!(fig)

    save(joinpath(supplementalspath, "hydro_$(lowercase(string(ssp))).png"), fig, px_per_unit = 3)
end

### Overlap of fasciola risk and livestock
livestocklazy = FD.get_livestock_data(; lazy = true)
livestock = disaggregate(read(crop(livestocklazy; to = ext, atol = 1e-4)), 2) # to match resolution of risk

fig6 = let fig = Figure(fontsize = 10)
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
                    FD.rasters_to_cmap(risk[ti][fs][Ti = 2, ssp = At(SSP370)], livestock[l_sp]; cmap)
                )           
                hidedecorations!(ax); hidespines!(ax)

                FD.Label_subplot(fig[l_idx, f_idx+ti_idx]; n = l_idx*4+f_idx+ti*2-6, valign = 1.02)

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
        ylabel = "Livestock density (/km²)", xlabel = rich(rich("Fasciola", font = :italic), " risk"),
        xticklabelrotation = pi/4,
        alignmode = Mixed(bottom = 0) # hack to include protrusions so we don't get excess white space
    )

    colsize!(fig.layout, 5, Relative(0.15))

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    resize_to_layout!(fig)
    fig 
end;
save(joinpath(figurepath, "figure6.png"), fig6, px_per_unit = 10)


### Overlap of Fasciola gigantica and Fasciola hepatica
fig5 = let fig = Figure(fontsize = 12)

    xticks = [0, 0.25, 0.5, 0.75, 1]
    yticks = copy(xticks)
    brewer_seqseq2 = [
        colorant"#ffac36" colorant"black"
        colorant"#f3f3f3" colorant"#209ebe"
    ]
    cmap = FD.bivariate_colormap(xticks, yticks; colors = brewer_seqseq2)

    ax_current = Axis(fig[1, 1])
    plot!(
        ax_current,
        FD.rasters_to_cmap(risk.current.hepatica, risk.current.gigantica; cmap)
    )  
    hidedecorations!(ax_current); hidespines!(ax_current)
    Label(fig[0, 1], "Current", tellwidth = false, font = :bold)
    FD.Label_subplot(fig[1, 1]; n = 1)

    for d in eachindex(dates)
        ax = Axis(fig[1, d+1])
        plot!(
            ax,
            FD.rasters_to_cmap(
                risk.future.hepatica[Ti = d, ssp = At(SSP370)], 
                risk.future.gigantica[Ti = d, ssp = At(SSP370)];
                cmap
            )
        )  
        hidedecorations!(ax); hidespines!(ax)
        # date label
        d1, d2 = val(Rasters.span(dates))[d, :]
        Label(fig[0, d+1], "$(year(d1))-$(year(d2))"; tellwidth = false, font = :bold)
        FD.Label_subplot(fig[1, d+1]; n = d+1)
    end
    rowsize!(fig.layout, 1, Aspect(1, 1/asp))
    resize_to_layout!(fig)

    FD.bivariate_cmap_legend(
        fig[1, 4], cmap; 
        aspect = AxisAspect(1),
        ylabel = rich("Fasciola gigantica", font =:italic), xlabel = rich("Fasciola hepatica"; font = :italic),
        xticklabelrotation = pi/4,
        alignmode = Mixed(bottom = 0), # hack to include protrusions so we don't get excess white space
    )

    Label(fig[1,4], "Transmission risk", valign = 0.95, font = :bold)

    colsize!(fig.layout, 4, Relative(0.15))

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    resize_to_layout!(fig)
    fig 
end
save(joinpath(figurepath, "figure5.png"), fig5; px_per_unit)


#########################
# Supplementary figures #
#########################

#### occurrence data
figs1 = let
    samplingcoordinates = mapreduce(x -> values.(x), vcat, samplingbg)
    snailcoordinates = map(x -> x.geo, ocs)

    fig = Figure(size = (500, 700))
    ax = Axis(
        fig[1, 1], 
        limits = (-25, 53, -35, 73), 
        title = "Recorded occurrences of Lymnaeids"
    ) # probably use geoaxis instead
    hidedecorations!(ax); hidespines!(ax)
    poly!(ax, countries.geometry, color = :transparent, strokewidth = 0.7)

    scatter!(samplingcoordinates; color = :grey, label = "Other Lymnaeids", markersize = 6)
    scatter!(snailcoordinates.radix; color = :red, label = rich("Radix natalensis", font = :italic), markersize = 6)
    scatter!(snailcoordinates.galba; color = :blue, label = rich("Galba truncatula", font = :italic), markersize = 6)
    axislegend(ax, position = (0.15, 0.15))
    colsize!(fig.layout, 1, Aspect(1, asp))
    resize_to_layout!(fig)

    fig
end;
save(joinpath(supplementalspath, "occurrences.png"), figs1; px_per_unit)

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


#### Predictions per model
for f in fasc_sps
    rfh = risk_all.future[f][Ti = 2, ssp = 2]
    fig_allmodels = let fig = Figure(fontsize = 8)
        for (i, gcm) in enumerate(dims(rfh, :gcm))
            for (j, ghm) in enumerate(dims(rfh, :ghm))
                ax = Axis(fig[j,i])
                plot!(ax, rfh[gcm = At(gcm), ghm = At(ghm)], colorrange = (0,1), colormap = Reverse(:Spectral))
                hidespines!(ax); hidedecorations!(ax)
                if i == 1
                    Label(fig[j,0], string(ghm), rotation = pi/2, tellheight = false)
                end
            end
            Label(fig[0,i], replace(string(gcm), "_" => "-"), tellwidth = false)
            colsize!(fig.layout, i, Aspect(1, asp))
        end
        colgap!(fig.layout, 3)
        rowgap!(fig.layout, 3)
        resize_to_layout!(fig)
        fig
        end
    save(joinpath(supplementalspath, "$(f)_by_gcm_ghm_SSP370_2081.png"), fig_allmodels)
end

# useful functions
function scatter_by_group!(ax, data_x, data_y, groups; marker = fill(:circle, length(data_x)), kwargs...)
    n_groups = length(unique(groups))

    for (i, group) in enumerate(unique(groups))
        indices = findall(groups .== group)
        scatter!(ax, data_x[indices], data_y[indices]; color = i, colorrange = (1, n_groups), label = group, marker = marker[indices], kwargs...)
    end
    return nothing
end

function plot_quantiles!(ax, temperatures, preds;  quantiles = [0.5, 0.025, 0.975], linestyle = [:solid, :dash, :dash], color = :black, kw...)
    for (q, linestyle) in zip(quantiles, linestyle)
        pred_quantiles = quantile.(eachcol(preds), q)
        lines!(ax, temperatures, pred_quantiles; linestyle, color, kw...)
    end
end

function plot_life_history(gl, data, temperatures; title, colormap = :rainbow, legend = true, kw...)
    units = unique(data.unit)
    xlabel = "temperature (°C)"

    if first(data.trait) == "infection efficiency" && first(data.fasciola) == "hepatica"
        # in this case, we use both info from the share of infected snails and radioactivity-based measurements
        # this means we need a double axis and a more complex legend
        data_N = filter(row -> row.unit == "N", data)
        shares = data_N.observation ./ data_N.sample_size
        data_rad = filter(row -> startswith(row.unit, "radioactivity"), data)

        @assert nrow(data) == nrow(data_N) + nrow(data_rad) # this should catch anything unexpected

        ax = Axis(gl[1, 1], limits = (5, 40, 0, 1); title, xlabel, ylabel = "share", aspect = 1, kw...)
        ax2 = Axis(gl[1, 1], limits = (5, 40, 0, 300); ylabel = "radioactivity (CPM)", 
            yaxisposition = :right, ygridvisible = false, aspect = 1, kw...)
        hidespines!(ax2); hidexdecorations!(ax2)

        scatter_by_group!(ax, data_N.temperature, shares, data_N.Citation; colormap)
        scatter_by_group!(ax2, data_rad.temperature, data_rad.observation, data_rad.Citation; colormap = :blues)
        
        if legend
            pl1 = Makie.get_labeled_plots(ax, merge = false, unique = true)
            pl2 = Makie.get_labeled_plots(ax2, merge = false, unique = true)
            Legend(gl[1,2], [pl1[1], pl2[1]], [pl1[2], pl2[2]], ["share", "radioactivity"], labelsize = 12)
        end

    else
        observations, ylabel = if units == ["N"]
            (data.observation ./ data.sample_size, "share")
        else
            (data.observation, units[1])
        end

        marker = if units == ["days"]
            ifelse.(.~ismissing.(data.exp_length) .&& data.exp_length .<= data.observation, :utriangle, :circle)
        else
            fill(:circle, nrow(data))
        end
        yscale = units == ["cercariae/snail"] ? log10 : identity
        ymin = yscale == identity ? 0 : 1
        ax = Axis(gl[1, 1], limits = (minimum(temperatures), maximum(temperatures), ymin, nothing); title, xlabel, ylabel, yscale, aspect = 1, kw...)
        scatter_by_group!(ax, data.temperature, observations, data.Citation; colormap, marker)
        if legend
            Legend(gl[1,2], ax;  labelsize = 12, orientation = :vertical)
        end
    end

    return ax

end
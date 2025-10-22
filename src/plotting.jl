# useful functions
function scatter_by_group!(ax, data_x, data_y, groups, labels=unique(groups); marker = fill(:circle, length(data_x)), kwargs...)
    n_groups = length(unique(groups))

    for (i, group) in enumerate(unique(groups))
        indices = findall(groups .== group)
        scatter!(
            ax, data_x[indices], data_y[indices]; 
            color = i, colorrange = (1, n_groups), 
            label = labels[i], marker = marker[indices], kwargs...
        )
    end
    return nothing
end

function plot_quantiles!(
    ax, temperatures, preds; 
    quantiles = [0.025, 0.975], 
    linestyle = :dash, 
    color = :black, 
    plot_mean = true,
    kw...)
    ls = map(quantiles) do q
        pred_quantiles = StatsBase.quantile.(eachcol(preds), q)
        lines!(ax, temperatures, pred_quantiles; linestyle, color, kw...)
    end
    if plot_mean
        mean_pred = vec(Statistics.mean(preds, dims = 1))
        ls_m = lines!(ax, temperatures, mean_pred; color, linestyle = :solid, kw...)
        ls = vcat(ls, ls_m)
    end
    return ls
end
function shade_quantiles!(ax, temperatures, preds; p, color, maxval = 1e5, minval = 1e-5)
    pred_quantiles = StatsBase.quantile.(eachcol(preds), Ref(p))
    upper = min.(last.(pred_quantiles), maxval)
    lower = min.(max.(first.(pred_quantiles), minval), upper)
    band!(
        ax, temperatures, 
        lower,
        upper;
        color
    )
end

function plot_life_history!(ax, data, temperatures; colormap = :rainbow, legend_gp = nothing)
    units = unique(data.unit)
    ax.xlabel = "temperature (°C)"
    ax.aspect = 1

    if first(data.trait) == "infection efficiency" && first(data.fasciola) == "hepatica"
        # in this case, we use both info from the share of infected snails and radioactivity-based measurements
        # this means we need a double axis and a more complex legend
        data_N = filter(row -> row.unit == "N", data)
        shares = data_N.observation ./ data_N.sample_size
        data_rad = filter(row -> startswith(row.unit, "radioactivity"), data)

        @assert nrow(data) == nrow(data_N) + nrow(data_rad) # this should catch anything unexpected

        ax.limits = (extrema(temperatures)..., 0, 1)
        ax.ylabel = "Share"

        # create a secondary axis for plotting the radioactivity data
        ax2 = Axis(gridpos(ax);
            limits = (extrema(temperatures)..., 0, nothing),
            xlabel = "temperature (°C)",
            ygridvisible = false,
            yaxisposition = :right,
            alignmode = ax.alignmode[],
            aspect = 1
        )
        hidespines!(ax2)#; hidexdecorations!(ax2)

        scatter_by_group!(ax, data_N.temperature, shares, data_N.Citation; colormap)
        scatter_by_group!(
            ax2, data_rad.temperature, data_rad.observation, 
            data_rad.group, data_rad.Citation; 
            colormap = :blues
        )
        
        if !isnothing(legend_gp)
            pl1 = Makie.get_labeled_plots(ax, merge = false, unique = false)
            pl2 = Makie.get_labeled_plots(ax2, merge = false, unique = false)
            Legend(
                legend_gp, [pl1[1], pl2[1]], [pl1[2], pl2[2]], ["share", "radioactivity"]; 
                labelsize = 12, tellheight = true)
        end
        return ax2

    else
        isshare = units == ["N"]
        observations, ax.ylabel = if isshare
            # convert absolute numbers to share for plotting
            (data.observation ./ data.sample_size, "share")
        else
            (data.observation, units[1])
        end

        marker = if units == ["days"]
            ifelse.(.~ismissing.(data.exp_length) .&& data.exp_length .<= data.observation, :utriangle, :circle)
        else
            fill(:circle, nrow(data))
        end
        if units == ["cercariae/snail"] 
            ymin = 1
            ymax = maximum(skipmissing(observations)) * 1.2
            ax.yscale = log10
        else
            ymin = 0
            ymax = units == ["N"] ? 1 : maximum(skipmissing(observations)) * 1.2
        end
        ax.limits = (extrema(temperatures)..., ymin, ymax)

        scatter_by_group!(ax, data.temperature, observations, data.Citation; colormap, marker)
        if !isnothing(legend_gp)
            Legend(legend_gp, ax; labelsize = 12, orientation = :vertical, tellheight = true)
        end
    end

    return ax

end
function gridpos(ax::Axis)
    sp = Makie.GridLayoutBase.gridcontent(ax).span
    si = Makie.GridLayoutBase.gridcontent(ax).side
    gl = Makie.GridLayoutBase.gridcontent(ax).parent
    GridPosition(gl, sp, si)
end


# bivariate color map utils
using Makie: Colors
import Makie.Colors: RGBA, RGB
function rasters_to_cmap(x, y; cmap, nacol = RGB(1, 1, 1))
    fallback = convert(eltype(cmap), nacol)
    function f(xval, yval)
        if ismissing(xval) || ismissing(yval) || xval === missingval(x) || yval === missingval(y) 
            fallback
        else
            cmap[X = Contains(xval), Y = Contains(yval)]
        end
    end
    # if only broadcast_dims had a strict=false
    @d f.(x, y) strict = false
end

function bivariate_cmap_legend(
    gp, cmap; 
    halign = :center, valign = :center, 
    height = nothing, width = nothing, 
    xticksf = string, yticksf = string,
    kw... # to pass to axis
)
    xspan, yspan = val.(Rasters.span(dims(cmap)))
    xticks = [xspan[1]; xspan[2,:]]
    yticks = [yspan[1]; yspan[2,:]]
    ax = Axis(
        gp; 
        halign, valign, height, width, 
        tellheight = false, tellwidth = false,
        xticks = (0:size(xspan,2), xticksf.(xticks)),
        yticks = (0:size(yspan,2), yticksf.(yticks)),
        kw...
    )
    heatmap!(ax, 0:size(xspan,2), 0:size(yspan,2), parent(cmap))
end

function bivariate_colormap(xticks, yticks; colors)
    palette = colors2d.(range(0,1, length = length(xticks)-1), collect(range(0,1, length = length(yticks)-1))'; colors)
    ds = map((xticks, yticks), (X,Y)) do ticks, D
        D(Sampled(ticks[1:end-1]; span = Explicit(rotl90([ticks[2:end] ticks[1:end-1]]))))
    end
    DimArray(palette, ds)
end

function colors2d(x, y; colors)
    x1 = colors[2,2] * x + colors[2,1] * (1-x)
    x2 = colors[1,2] * x + colors[1,1] * (1-x)
    return x2 * y + x1 * (1-y)
end

function Label_subplot(gp; n, kw...)
    text = "$('A'+(n-1)))"
    Label(
        gp, text, font = :bold, tellheight = false, tellwidth = false, 
        halign = :left, valign = 1.02; kw...
    )
end

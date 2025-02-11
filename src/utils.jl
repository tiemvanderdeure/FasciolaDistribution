# Just a convenient function because namedtuples of namedtuples are used throughout
mapmap(f, x) = map(y -> map(f, y), x)

function quarterly_means(transmission_raster)
    grouped = Rasters.groupby(transmission_raster, Dim{:month} => CyclicBins(identity;cycle = 12, start = 12, step = 3))
    transmission_raster_quarterly = dropdims.(mean.(grouped; dims = Dim{:month}); dims = Dim{:month})
    transmission_raster_quarterly = cat(transmission_raster_quarterly...; dims = Dim{:quarter}(1:4))
end

# to be able to write rasters
function writeable_dims(A)
    rebuild(A; dims = writeable_dims(dims(A)))
end

function writeable_dims(ds::Rasters.DD.DimTuple)
    for s in (:ssp, :gcm)
        d = dims(ds, s)
        ds = isnothing(d) ? ds : set(ds, d => string.(d))
    end
    d = dims(ds, Ti)
    ds = if isnothing(d) 
        ds 
    else
        # cannot write dates so convert to julian - need to convert both span and actual lookups 
        ds = set(ds, d => Explicit(Dates.datetime2julian.(DateTime.(val(Lookups.span(d))))))
        set(ds, d => Dates.datetime2julian.(DateTime.(d)))
    end
    return ds
end

function reformat_dims(A::AbstractRaster)
    # format these back from string to datatypes
    for s in (:ssp, :gcm)
        d = dims(A, s)
        A = isnothing(d) ? A : set(A, d => eval.(Symbol.(d)))
    end
    # format time back from julian to datetime
    d = dims(A, Ti)
    A = if isnothing(d) 
        A 
    else
        set(A, d => Date.(Dates.julian2datetime.(d)))
        set(A, d => Explicit(Date.(Dates.julian2datetime.(val(Lookups.span(d))))))
    end
    return A
end

# See https://github.com/JuliaGeo/NCDatasets.jl/issues/274
Base.delete_method.(methods(view, (NCDatasets.CommonDataModel.AbstractVariable, Colon)))

# To make lazy arrays work with @d, see https://github.com/rafaqz/DimensionalData.jl/issues/925
# This is all copy-pasted from the @d macro except for the last bit.
const DD = DimensionalData
macro lazyd(expr::Expr, options::Union{Expr,Nothing}=nothing)
    options_dict, options_expr = DD._process_d_macro_options(options)
    broadcast_expr, var_list = DD._find_broadcast_vars(expr)
    var_list_assignments = map(var_list) do (name, expr)
        Expr(:(=), name, expr)
    end
    vars_expr = esc(Expr(:tuple, map(first, var_list)...))
    var_list_expr = esc(Expr(:block, var_list_assignments...))
    dims_expr = if haskey(options_dict, :dims)
        order_dims = options_dict[:dims]
        quote
            order_dims = $order_dims
            found_dims = DD._find_dims(vars)
            all(hasdim(order_dims, found_dims)) || 
                throw(ArgumentError("order $(DD.basedims(order_dims)) dont match dimensions found in arrays $(DD.basedims(found_dims))"))
            dims = $DimensionalData.dims(found_dims, order_dims)
        end
    else
        quote
            dims = DD._find_dims(vars)
        end
    end
    # only this is different from @d
    rebuild_expr = quote
        bc = Base.Broadcast.instantiate(DD._unwrap_broadcasted(LazyArrays.lazy.($broadcast_expr)))
        rebuild(first(vars); data = LazyArray(bc), dims)
    end
    quote
        let
            options = $options_expr
            $var_list_expr
            vars = $vars_expr
            $dims_expr
            $rebuild_expr
        end
    end
end

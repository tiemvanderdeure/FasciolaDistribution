using Rasters.Lookups: Categorical
# SSPs to use
ssps = Categorical([SSP126, SSP370]) |> Dim{:ssp}
# Global climate models
gcms = Categorical([GFDL_ESM4, IPSL_CM6A_LR, MPI_ESM1_2_HR, UKESM1_0_LL, MRI_ESM2_0]) |> Dim{:gcm} # 
# Global hydrological models
ghms = Categorical(["CWatM", "JULES-W2", "MIROC-INTEG-LAND", "WaterGAP2-2e"]) |> Dim{:ghm} # 
# timepoints
datematrix = [
    Date(2041) Date(2061) - Day(1)
    Date(2081) Date(2101) - Day(1)
]
dates = Ti(Sampled(datematrix[:,1]; span=Explicit(datematrix))) |> Rasters.DD.format


function get_countries(continents)
    countries = naturalearth("admin_0_countries", 10)
    countries.geometry[in.(countries.CONTINENT, Ref(continents))]
end

function crop_and_read(ras, geoms)
    ras_crop = read(crop(ras; to = geoms))
   # ras_roi = mask(ras_crop; with = geoms, boundary = :touches)
   # Rasters.trim(ras_roi)
end

function get_temperature_data(geoms; gcms, ssps, dates)
    tas = Rasters.RasterSeries(WorldClim{Climate}, :tavg; month = 1:12, res = "2.5m", lazy = true, missingval = missingval) |> Rasters.combine
    tas = crop_and_read(tas, geoms)
    tminmaxf = (@d _get_future_temp.(gcms, ssps, dates)) |> RasterSeries |> Rasters.combine
    tminmaxf = crop_and_read(tminmaxf, geoms)
    tasf = middle.(tminmaxf.tmin, tminmaxf.tmax)
    tasf = set(tasf, Rasters.Band => Rasters.format(Dim{:month}(1:12)))
    return (current = tas, future = tasf)
end
function get_bioclim(layers::NTuple{<:Any,Int}, geoms; gcms, ssps, dates)
    bio = RasterStack(WorldClim{BioClim}, layers; res = "2.5m", lazy = true)
    biof = (@d _get_future_worldclim.(gcms, ssps, dates)) |> RasterSeries |> Rasters.combine
    biof = view(biof, Band = collect(layers))
    bio = crop_and_read(bio, geoms)
    biof = crop_and_read(biof, geoms)
    biof = RasterStack(biof; layersfrom = Rasters.Band, name = keys(bio))
    return (current = bio, future = biof)
end

_get_future_temp(gcm, ssp, date) =
    RasterStack(WorldClim{Future{Climate, CMIP6, gcm, ssp}}, (:tmin,:tmax); res = "2.5m", date, lazy = true, missingval = missingval)
_get_future_worldclim(gcm, ssp, date) =
    Raster(getraster(WorldClim{Future{BioClim, CMIP6, gcm, ssp}}; res = "2.5m", date), lazy = true)



# current, only ghm
function discharge_url(ghm)
    URIs.URI(
        scheme = "https",
        host = "files.isimip.org",
        path = "/ISIMIP3a/OutputData/water_global/$ghm/gswp3-w5e5/historical/" * discharge_filename(ghm)
    )
end
# future, with gcm and ssp
function discharge_url(gcm, ssp, ghm)
    _gcm = lowercase(RDS._format(gcm))
    URIs.URI(
        scheme = "https",
        host = "files.isimip.org",
        path = "/ISIMIP3b/OutputData/water_global/$ghm/$_gcm/future/" * discharge_filename(gcm, ssp, ghm)
    )
end
function discharge_filepath(ghm)
    joinpath(
        RDS.rasterpath(),
        "ISIMIP",
        "water_global",
        "Historical",
        ghm,
        discharge_filename(ghm)
    )
end
function discharge_filepath(gcm, ssp, ghm)  
    joinpath(
        RDS.rasterpath(), 
        "ISIMIP", 
        "water_global", 
        "Future", 
        RDS._format(ssp), 
        RDS._format(gcm), 
        ghm, 
        discharge_filename(gcm, ssp, ghm)
    )
end
# current
discharge_filename(ghm) = "$(lowercase(ghm))_gswp3-w5e5_obsclim_histsoc_default_dis_global_monthly_1901_2019.nc"
function discharge_filename(gcm, ssp, ghm)
    _ssp = lowercase(RDS._format(ssp))
    _gcm = lowercase(RDS._format(gcm))
    _ghm = lowercase(ghm)
    "$(_ghm)_$(_gcm)_w5e5_$(_ssp)_2015soc-from-histsoc_default_dis_global_monthly_2015_2100.nc"
end
function discharge_rasterpath(args...)
    url = discharge_url(args...)
    filepath = discharge_filepath(args...)
    RDS._maybe_download(url, filepath)
    return filepath
end

function get_discharge_data(geoms; gcms, ssps, ghms, dates)
    discharge_c = Rasters.RasterSeries((@d discharge_rasterpath.(ghms)); lazy = true) |> Rasters.combine
    discharge_f = Rasters.RasterSeries((@d discharge_rasterpath.(gcms, ssps, ghms)); lazy = true) |> Rasters.combine
    # fix the dimensions
    current, future = map((discharge_c, discharge_f)) do d
        d = set(d, 
            X => Intervals(Center()), 
            Y => Intervals(Center()), 
            Ti => Sampled(sampling = Intervals(Start()), span = Regular(Month(1)))
        )
    end
    current = crop_and_read(view(current, Ti(Date(1980) .. Date(2010))), geoms)

    monthly_current = groupbymonth(current)
    monthly_future = groupbymonth.(crop_and_read.(groupby(future, dates), Ref(geoms)))
    monthly_future = Rasters.combine(RasterSeries(monthly_future))

    return (; current = monthly_current, future = monthly_future)
end

function groupbymonth(x)
    groups = groupby(x, Ti => month)

    monthly_discharge = map(groups) do g
        # this super convoluted thing is necessary because some of the rasters are filled with missing values
        # in particular the last time slice of the JULES mdoel
        mean.(skipmissing.(eachslice(g, dims = otherdims(x, Ti))))
    end |> RasterSeries |> Rasters.combine
    # Time dimension is now month
    monthly_discharge = set(monthly_discharge, Ti => Dim{:month}(1:12))
    # missings became NaNs
    replace_missing(monthly_discharge, NaN)
end

using NearestNeighbors, StaticArrays
function water_distance(ras)
    rivers = naturalearth("ne_10m_rivers_lake_centerlines")
    rivergeoms =filter(!isempty , rivers.geometry)
    lakes = naturalearth("ne_10m_lakes")
    lakegeoms = lakes.geometry
    bm = boolmask(ras)
    riverraster = rasterize(last, rivergeoms; to = bm, missingval = false, fill = true)
    lakesraster = rasterize(last, lakes; to = bm, missingval = false, fill = true)
    waterraster = riverraster .|| lakesraster
    points = Rasters.DimPoints(waterraster) .|> SVector
    waterpoints = points[waterraster]

    water_dist_rast = Raster{Union{Missing, Float32}}(undef, dims(ras, (X,Y)); name = :water_dist)
    fill!(water_dist_rast, missing)

    kdtree = KDTree(waterpoints; leafsize = 20)

    let kdtree = kdtree
        map!(view(water_dist_rast, bm), points[bm]) do point
            knn(kdtree, point, 1)[2][1]
        end
    end

    return water_dist_rast
end
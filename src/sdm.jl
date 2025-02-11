const GBIFids_files = Dict(
    "lymnaeids_europe" => "0000023-241007104925546",
    "lymnaeids_africa" => "0000137-241007104925546"
)

function get_snail_occurrences()
    _maybe_download_gbif()

    museum = CSV.read("data/snails/pulmonates.csv", DataFrame; delim = ";", types = Dict(18 => String))
    museum.Longitude .= map(r -> ismissing(r) || contains("°")(r) ? missing : r, museum.Longitude)
    museum.Latitude .= map(r -> ismissing(r) || contains("°")(r) ? missing : r, museum.Latitude)
    dropmissing!(museum, [:Latitude, :Longitude, :Species, :Genus])
    museum.X = parse.(Float64, museum.Longitude)
    museum.Y = parse.(Float64, museum.Latitude)
    radix_museum = filter(m -> m.Species == "natalensis" && contains("Lymnaea")(m.Genus), museum)
    galba_museum = filter(row -> row.Species == "truncatula", museum)

    lymnaeids_europe = process_gbif("data/snails/lymnaeids_europe.csv")
    lymnaeids_africa = process_gbif("data/snails/lymnaeids_africa.csv")

    # We do some harsher filtering for europe
    filter!(row -> row.year > 1980 && (ismissing(row.coordinateUncertaintyInMeters) || row.coordinateUncertaintyInMeters < 1000), lymnaeids_europe)
    # almost all records from ireland are from before 1970 and have coordinateUncertaintyInMeters ≈ 7000
    # For africa we have so few G. trunctula records that we just keep them all

    galba_gbif_af = filter(row -> row.species == "Galba truncatula", lymnaeids_africa)
    galba_gbif_eu = filter(row -> row.species == "Galba truncatula", lymnaeids_europe)
    radix_gbif = filter(row -> row.species == "Radix natalensis", lymnaeids_africa)

    radix = [xy_from_df(radix_museum); xy_from_df(radix_gbif)]
    galba_af = [xy_from_df(galba_museum); xy_from_df(galba_gbif_af)]
    galba_eu = xy_from_df(galba_gbif_eu)

    sampling_africa = [xy_from_df(museum); xy_from_df(lymnaeids_africa)]
    sampling_eu = xy_from_df(lymnaeids_europe)

    sampling_radix = setdiff(sampling_africa, radix)
    sampling_galba_eu = setdiff(sampling_eu, galba_eu)
    sampling_galba_af = setdiff(sampling_africa, galba_af)

    return ((; galba_eu, galba_af, radix), (; galba_eu = sampling_galba_eu, galba_af = sampling_galba_af, radix = sampling_radix))
end

function background_sample(r, i)
    bm = Rasters.boolmask(r)
    sample(r, weight = weights(bm), i)
end

function process_gbif(file)
    df = CSV.read(file, DataFrame; delim ='\t', select = [:decimalLatitude, :decimalLongitude, :year, :species, :coordinateUncertaintyInMeters, :countryCode])
    rename!(df, :decimalLatitude => :Y, :decimalLongitude => :X)
    dropmissing!(df, [:X, :Y, :year, :species])
    return df
end

xy_from_df(df) = (X = df.X, Y = df.Y) |> Tables.rowtable |> unique

function mypredict(m, bio)
    bm = Rasters.boolmask(bio)
    out = similar(first(layers(bio)))
    out = rebuild(out; missingval = eltype(out)(NaN))
    fill!(out, NaN)
    @views out[bm] .= predict(m, bio[bm])
    return out
end

function _maybe_download_gbif()
    snaildatadir = joinpath("data", "snails")
    for file in keys(GBIFids_files)
        filenamebase = joinpath(snaildatadir, file)
        if !isfile(filenamebase * ".csv")
            GBIF2.occurrence_download(GBIFids_files[file]; filename = filenamebase * ".zip")
            write(filenamebase * ".csv", read(ZipFile.Reader(filenamebase * ".zip").files[1]))
        end
    end
end
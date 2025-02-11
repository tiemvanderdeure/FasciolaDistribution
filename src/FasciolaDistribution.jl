__precompile__(false)
module FasciolaDistribution

using NaturalEarth, Rasters, RasterDataSources, ArchGDAL, NCDatasets, Rasters.Lookups
using Dates, URIs, ZipFile, LazyArrays
using CSV, DataFrames, GBIF2
using GLM, Turing, StatsBase, Random
using Makie, BibParser

const RDS = RasterDataSources

const traits_str = ["hatching time", "prepatent period", "infection efficiency", "cercarial release", "hatching success"]
const trait_keys = Tuple(Symbol.(replace.(traits_str, " " => "_")))

export traits_str, trait_keys, load_life_history_data, define_life_history_models, life_cycle_model, find_max_life_cycle
export scatter_by_group!, plot_quantiles!, plot_life_history
export get_temperature_data, get_future_temperature_data, get_discharge_data, get_bioclim
export get_snail_occurrences
export quarterly_means, mapmap, writeable_dims, lazyd
export mypredict


include("utils.jl")
include("models.jl")
include("fits.jl")
include("life_cycle_model.jl")
include("plotting.jl")
include("climate.jl")
include("sdm.jl")

end
using SummaryTables, WriteDocx, Printf
global supplementary_tablespath = joinpath("images", "supplementary")

### Table with studies and what life history traits they contribute with

# get a vector of keys for each trait-fasciola species combination
trait_to_key = mapmap(life_history_data) do data
    unique(data.Key)
end

# get the total of all keys 
allkeys = unique(vcat(trait_to_key.hepatica..., trait_to_key.gigantica...))
keys_to_traits = Dict([k => String[] for k in allkeys])
keys_to_species = Dict([k => String[] for k in allkeys])

key_to_trait = map(keys(trait_to_key)) do f_sp
    map(keys(trait_to_key[f_sp])) do trait
        for key in trait_to_key[f_sp][trait]
            push!(keys_to_traits[key], string(trait))
            push!(keys_to_species[key], string(f_sp))
        end
    end
end

# combine into a single Dataframe
df_sp = DataFrame(Key = collect(keys(keys_to_species)), Species = unique!.(collect(values(keys_to_species))))
df_traits = DataFrame(Key = collect(keys(keys_to_traits)), Trait = unique!.(collect(values(keys_to_traits))))
df_citations = FasciolaDistribution.life_history_citations()

df = leftjoin(leftjoin(df_sp, df_traits, on = :Key), df_citations, on = :Key)
# format species and trait
df.Species = join.(df.Species, " & ")
df.Trait = join.(df.Trait, ", ")

# conver to a Table
cells1 = hcat(Cell.(df.Citation), Cell.(df.Species), Cell.(df.Trait))
table1 = Table([Cell.(["Citation", "Species", "Trait"])'; cells1], header = 1)

### Table with posterior distributions for supplementary materials

# helper functions
argnames(fix::FasciolaDistribution.FixN) = 
    Base.method_argnames(first(methods(fix.f, (map(typeof, fix.args)..., Float64))))[2:end-1]
function format_posterior(d::Vector{<:Real}) 
    q_025, q_975 = quantile(d, [0.025, 0.975])
    m = mean(d)
    nd = max(2 - floor(Int, log(10, m)), 1)
    @sprintf "%.*f (%.*f-%.*f)" nd m nd q_025 nd q_975
end

format_prior(d) = replace(string(d), r"\{(.*)\}" => "")

posterior_tables = map(trait_keys) do trait
    @show trait
    names = argnames(first(first(gqs)[trait]))
    param_cells = Cell.(string.(names))
    # get a matrix of cells - for hepatica and gigantica
    datacells = mapreduce(hcat, (:hepatica, :gigantica)) do fasciola
        posteriors_raw = getfield.(gqs[fasciola][trait], :args)
        posteriors = NamedTuple(K => getindex.(posteriors_raw, i) for (i, K) in enumerate(names))
        collect(map(p -> Cell(format_posterior(p)), posteriors))
    end

    priors = life_history_models.hepatica[trait].defaults.priors
    priorcells = [Cell(format_prior(priors[n])) for n in names]

    top_row = Cell[
        Cell("Parameter"; bold = true)
        Cell("F. hepatica"; bold=true, italic=true)
        Cell("F. gigantica"; bold=true, italic=true)
        Cell("Prior"; bold=true)
    ] |> permutedims

    Table([top_row; [param_cells datacells priorcells]])
end |> NamedTuple{trait_keys}

tables = mapreduce(vcat, enumerate(trait_keys)) do (i, K)
    [
        WriteDocx.Paragraph([
            WriteDocx.Run([WriteDocx.Text("Table $i: $K")]),
        ]),
        SummaryTables.to_docx(posterior_tables[K])
    ]
end |> WriteDocx.Section

tables = WriteDocx.Section(SummaryTables.to_docx.(collect(posterior_tables)))

WriteDocx.save(
    joinpath(supplementary_tablespath, "posterior_tables.docx"), 
    WriteDocx.Document(WriteDocx.Body([tables]))
)

supplementals = W.Document(
    W.Body([
        W.Section([
            W.Paragraph([
                W.Run([W.Text("Table 1")]),
            ]),
            SummaryTables.to_docx(table1),
            W.Paragraph([
                W.Run([W.Text("Table 2")]),
            ]),
            SummaryTables.to_docx(table2),
        ]),
    ]),
)

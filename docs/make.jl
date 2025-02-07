push!(LOAD_PATH, "../src/")
using Documenter, Agate

using Literate

using OceanBioME
using OceanBioME: SugarKelp, LOBSTER, NutrientPhytoplanktonZooplanktonDetritus
using OceanBioME.Sediments: SimpleMultiG, InstantRemineralisation
using OceanBioME: CarbonChemistry, GasExchange

using Oceananigans.Grids: RectilinearGrid

using CairoMakie
CairoMakie.activate!(type = "svg")

# Examples

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

examples = [
    "Size structures" => "size_structure",
]

example_scripts = [ filename * ".jl" for (title, filename) in examples ]

replace_silly_warning(content) = replace(content, r"┌ Warning:.*\s+└ @ JLD2 ~/\.julia/packages/JLD2/.*/reconstructing_datatypes\.jl.*\n" => "")

for example in example_scripts
    example_filepath = joinpath(EXAMPLES_DIR, example)

    withenv("JULIA_DEBUG" => "Literate") do
        Literate.markdown(example_filepath, OUTPUT_DIR; 
                          flavor = Literate.DocumenterFlavor(),
                          repo_root_url = "https://agate-model.github.io/Agate.jl/",
                          execute = true,
                          postprocess = replace_silly_warning)
    end
end

example_pages = [ title => "generated/$(filename).md" for (title, filename) in examples ]


makedocs(;
    sitename="Agate.jl Documentation",
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    # modules=[Agate],
    pages=[
        "About" => "index.md",
        "API" => "api.md",
        "Library" => "library.md",
        "Examples" => example_pages,
    ],
)

deploydocs(; repo="https://github.com/agate-model/Agate.jl.git")

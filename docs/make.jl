push!(LOAD_PATH, "../src/")
using Documenter, Literate
using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using CairoMakie

# Examples

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR = joinpath(@__DIR__, "src/generated")

examples = [
    "Box model factories" => "box_model_factories",
    "Box model: add heterotroph (H)" => "box_model_add_heterotroph",
    "Column model" => "1D_column",
]

example_scripts = [filename * ".jl" for (title, filename) in examples]

function replace_silly_warning(content)
    return replace(
        content,
        r"┌ Warning:.*\s+└ @ JLD2 ~/\.julia/packages/JLD2/.*/reconstructing_datatypes\.jl.*\n" => "",
    )
end

for example in example_scripts
    example_filepath = joinpath(EXAMPLES_DIR, example)

    withenv("JULIA_DEBUG" => "Literate") do
        Literate.markdown(
            example_filepath,
            OUTPUT_DIR;
            flavor=Literate.DocumenterFlavor(),
            repo_root_url="https://oceanbiome.github.io/OceanBioME.jl",
            execute=true,
            postprocess=replace_silly_warning,
        )
    end
end

example_pages = [title => "generated/$(filename).md" for (title, filename) in examples]

model_pages = ["NiPiZD" => "nipizd.md", "DARWIN" => "darwin.md"]

makedocs(;
    sitename="Agate.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        size_threshold=500000,  # 500KB threshold
        size_threshold_warn=200000,
    ), # 200KB warning
    # modules=[Agate],
    pages=[
        "About" => "index.md",
        "Quick start" => "quick_start.md",
        "Models" => model_pages,
        "Library" => "library.md",
        "Examples" => example_pages,
        "API" => ["Constructor API" => "api_constructor.md", "Equation-based dynamics" => "equations_api.md", "Reference" => "api.md"],
    ],
)

deploydocs(; repo="https://github.com/agate-model/Agate.jl.git")

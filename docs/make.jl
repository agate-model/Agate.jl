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
    "Box model customisation" => "box_model_factories", "Column model" => "1D_column"
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
            repo_root_url="https://github.com/agate-model/Agate.jl",
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
        "Home" => "index.md",
        "Beginner guide" => [
            "Start here" => "beginner.md",
            "Quick start" => "quick_start.md",
            "Using built-in models" => model_pages,
            "Examples" => example_pages,
        ],
        "Developer guide" => [
            "Start here" => "developer.md",
            "Variants" => "variants.md",
            "Adding a model" => "adding_models.md",
            "Constructor API" => "api_constructor.md",
            "Parameters & interactions" => "equations_api.md",
            "Callable dynamics" => "functors_api.md",
            "API reference" => "api.md",
        ],
        "Library" => "library.md",
    ],
)

deploydocs(; repo="https://github.com/agate-model/Agate.jl.git")

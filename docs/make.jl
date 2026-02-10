# to build docs from terminal `julia --project=docs docs/make.jl`
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
    "Column model" => "1D_column"
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
            # Avoid JuliaFormatter rewriting NamedTuple literals like `(; x=1)` into invalid syntax.
            format=false,
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
    modules=[Agate],
    checkdocs=:exports,
    warnonly=[:missing_docs],
    pages=[
        "Home" => "index.md",
        "Quick start" => "quick_start.md",
        "Examples" => example_pages,
        "Models" => model_pages,
        "Library" => "library.md",
        "Implementing new models" => [
            "Overview" => "developer.md",
            "Variants" => "variants.md",
            "Adding a model" => "adding_models.md",
            "Constructor API" => "api_constructor.md",
            "Parameters & interactions" => "parameters_and_interactions.md",
            "Callable dynamics" => "callable_dynamics.md",
        ],
        "API reference" => "api.md",
    ],
)

deploydocs(; repo="https://github.com/agate-model/Agate.jl.git")
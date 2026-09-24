# to build docs from terminal `julia --project=docs docs/make.jl`
import Pkg
Pkg.activate(@__DIR__)
Pkg.develop(Pkg.PackageSpec(; path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

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

const BUILD_COLUMN_EXAMPLE = lowercase(get(ENV, "AGATE_DOCS_BUILD_COLUMN", get(ENV, "CI", "false") == "true" ? "false" : "true")) in ("1", "true", "yes")

examples = [
    "Size structure" => "size_structure",
    "Predator-prey palatability" => "predator_prey_palatability",
    "Comparing phytoplankton light strategies" => "named_pfts",
    "Allometric parameters" => "allometric_relationships",
    "Exporting a model definition" => "export_model_recipe",
]

if BUILD_COLUMN_EXAMPLE
    push!(examples, "Column model" => "1D_column")
end

defining_models_examples = [
    "Model with mixotrophy" => "mixotrophy",
    "Model with bacterioplankton" => "detritus_bacteria",
]

differentiable_modelling = [
    "Forward-mode AD sensitivity" => "forward_mode_ad_nipizd_sensitivity",
    "Reverse-mode AD sensitivity" => "reverse_mode_ad_nipizd_sensitivity",
]

example_scripts = [
    filename * ".jl" for (title, filename) in
    vcat(examples, defining_models_examples, differentiable_modelling)
]

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
defining_models_pages = [
    title => "generated/$(filename).md" for (title, filename) in defining_models_examples
]
differentiable_modelling_pages = [
    title => "generated/$(filename).md" for (title, filename) in differentiable_modelling
]

contributor_pages = ["Architecture" => "architecture_overview.md"]
model_pages = ["NiPiZD" => "nipizd.md"]

makedocs(;
    sitename="Agate.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        size_threshold=1_000_000,
        size_threshold_warn=300_000,
    ),
    modules=[Agate],
    checkdocs=:exports,
    warnonly=[:missing_docs],
    pages=[
        "Home" => "index.md",
        "Quick start" => "quick_start.md",
        "Examples" => example_pages,
        "Models" => model_pages,
        "Defining new models" => defining_models_pages,
        "Library" => "library.md",
        "Differentiable modelling" => differentiable_modelling_pages,
        "Contributor guide" => contributor_pages,
        "API reference" => "api.md",
    ],
)

deploydocs(; repo="https://github.com/agate-model/Agate.jl.git")

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

examples = [
    "Size structure" => "size_structure",
    "Predator-prey palatability" => "predator_prey_palatability",
    "Traceable box model" => "traceability_box_model",
    "Column model" => "1D_column",
]
dev = ["Variants" => "variant"]
differentiable_modelling = ["Forward-mode AD sensitivity" => "forwarddiff_nipizd_ode"]

example_scripts = [
    filename * ".jl" for (title, filename) in vcat(examples, dev, differentiable_modelling)
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
dev_pages = [title => "generated/$(filename).md" for (title, filename) in dev]
differentiable_modelling_pages = [
    title => "generated/$(filename).md" for (title, filename) in differentiable_modelling
]

contributor_pages = ["Architecture" => "architecture_overview.md"]
model_pages = ["NiPiZD" => "nipizd.md", "DARWIN" => "darwin.md"]

makedocs(;
    sitename="Agate.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        size_threshold=1_000_000,  # allow larger example pages with embedded figures
        size_threshold_warn=300_000,
    ), # 300KB warning
    modules=[Agate],
    checkdocs=:exports,
    warnonly=[:missing_docs],
    pages=[
        "Home" => "index.md",
        "Quick start" => "quick_start.md",
        "Examples" => example_pages,
        "Models" => model_pages,
        "Library" => "library.md",
        "Implementing new models" => dev_pages,
        "Differentiable modelling" => differentiable_modelling_pages,
        "Contributor guide" => contributor_pages,
        "API reference" => "api.md",
    ],
)

deploydocs(; repo="https://github.com/agate-model/Agate.jl.git")

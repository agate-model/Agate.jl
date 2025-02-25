push!(LOAD_PATH, "../src/")
using Documenter, Agate
model_pages = ["NiPiZD" => "nipizd.md", "DARWIN" => "darwin.md"]

makedocs(;
    sitename="Agate.jl",
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    # modules=[Agate],
    pages=[
        "About" => "index.md",
        "Models" => model_pages,
        "Library" => "library.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ],
)

deploydocs(; repo="https://github.com/agate-model/Agate.jl.git")

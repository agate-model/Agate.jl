push!(LOAD_PATH, "../src/")
using Documenter, Agate

makedocs(;
    sitename="Agate.jl Documentation",
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    # modules=[Agate],
    pages=[
        "About" => "index.md",
        "API" => "api.md",
        "Library" => "library.md",
        "Examples" => "examples.md",
    ],
)

deploydocs(; repo="https://github.com/agate-model/Agate.jl.git")

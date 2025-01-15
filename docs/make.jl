push!(LOAD_PATH, "../src/")
using Documenter, Agate

# by default `source="src"`
makedocs(;
    sitename="Agate.jl Documentation",
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    # modules=[Agate],
    pages=[
        "About" => "index.md",
    ]
)

deploydocs(; repo="https://github.com/agate-model/Agate.jl.git")

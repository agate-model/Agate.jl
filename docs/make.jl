using Documenter, Agate

# by default `source="src"`
makedocs(;
    sitename="Agate.jl Documentation",
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    modules=[Agate],
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(; repo="https://github.com/agate-model/Agate.jl.git")

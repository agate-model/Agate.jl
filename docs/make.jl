push!(LOAD_PATH, "../src/")
using Documenter, Agate

About = "Introduction" => "index.md"

# DevDocs = ["Documentation of Agate.jl internals" => ["devdocs/growth.md"]]

PAGES = [About]

makedocs(;
    sitename="Agate.jl",
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
)

deploydocs(; repo="https://github.com/agate-model/Agate.jl.git")

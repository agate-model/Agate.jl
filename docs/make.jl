push!(LOAD_PATH,"../src/")
using Documenter, Agate

About = "Introduction" => "index.md"

DevDocs = [
	"Documentation of Agate.jl internals" => [
		"devdocs/growth.md"
	]
]

PAGES = [
    About,
    DevDocs
    ]

makedocs(
    sitename="Agate.jl Documentation",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "https://github.com/agate-model/Agate.jl.git",
)

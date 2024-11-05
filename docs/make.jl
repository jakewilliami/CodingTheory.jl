using CodingTheory, Documenter

Documenter.makedocs(
    clean = true,
    doctest = true,
    modules = Module[CodingTheory],
    repo = "",
    highlightsig = true,
    sitename = "CodingTheory Documentation",
    expandfirst = [],
    pages = [
        "Index" => "index.md",
    ]
)

deploydocs(;
    repo  =  "github.com/jakewilliami/CodingTheory.jl.git",
)

# deploydocs(
#     target = "build",
#     repo   = "github.com/jakewilliami/CodingTheory.jl.git",
#     branch = "gh-pages",
#     devbranch = "master",
#     devurl = "dev",
#     versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
#     push_preview    = false
# )

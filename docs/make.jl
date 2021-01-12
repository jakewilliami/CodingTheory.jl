include(joinpath(dirname(@__DIR__), "src", "CodingTheory.jl"))
using Documenter, .CodingTheory

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

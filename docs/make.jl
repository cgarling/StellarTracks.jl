using StellarTracks
using Documenter
using DocumenterCitations: CitationBibliography

DocMeta.setdocmeta!(StellarTracks, :DocTestSetup, :(using StellarTracks); recursive=true)

# The `format` below makes it so that urls are set to "pretty" if you are pushing them to a hosting service, and basic if you are just using them locally to make browsing easier.

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric) # style=:authoryear

makedocs(;
    modules=[StellarTracks],
    authors="Chris Garling",
    sitename="StellarTracks.jl",
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical="https://cgarling.github.io/StellarTracks.jl",
        edit_link="main",
        assets=String["assets/citations.css"],
    ),
    pages=[
        "index.md",
        "parsec.md",
        "api.md",
        "refs.md",
        "doc_index.md"
    ],
    pagesonly=true,
    doctest=false,
    linkcheck=true,
    warnonly=[:missing_docs, :linkcheck],
    plugins=[bib]
)

deploydocs(;
    repo="github.com/cgarling/StellarTracks.jl",
    versions = ["stable" => "v^", "v#.#"],
    devbranch="main",
    push_preview=true
)

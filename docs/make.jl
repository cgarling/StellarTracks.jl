using StellarTracks
using Documenter

DocMeta.setdocmeta!(StellarTracks, :DocTestSetup, :(using StellarTracks); recursive=true)

# The `format` below makes it so that urls are set to "pretty" if you are pushing them to a hosting service, and basic if you are just using them locally to make browsing easier.

makedocs(;
    modules=[StellarTracks],
    authors="Chris Garling",
    sitename="StellarTracks.jl",
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical="https://cgarling.github.io/StellarTracks.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "PARSEC" => "parsec.md"
    ],
    doctest=false,
    linkcheck=true,
    warnonly=[:missing_docs, :linkcheck],
)

deploydocs(;
    repo="github.com/cgarling/StellarTracks.jl",
    versions = ["stable" => "v^", "v#.#"],
    devbranch="main",
    push_preview=true
)

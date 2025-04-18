using StellarTracks
using Documenter

DocMeta.setdocmeta!(StellarTracks, :DocTestSetup, :(using StellarTracks); recursive=true)

makedocs(;
    modules=[StellarTracks],
    authors="cgarling <chris.t.garling@gmail.com> and contributors",
    sitename="StellarTracks.jl",
    format=Documenter.HTML(;
        canonical="https://cgarling.github.io/StellarTracks.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    doctest=false,
    linkcheck=true,
    warnonly=[:missing_docs, :linkcheck],
)

deploydocs(;
    repo="github.com/cgarling/StellarTracks.jl",
    devbranch="main",
)

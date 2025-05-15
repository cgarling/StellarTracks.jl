using StellarTracks
using Documenter
using DocumenterCitations: CitationBibliography
using DocumenterInterLinks: InterLinks
using Unicode: normalize # for parsing CI environment variable

# Check if we are running on CI
if "CI" in keys(ENV)
    ci = parse(Bool, normalize(ENV["CI"]))
else
    ci = false
end

DocMeta.setdocmeta!(StellarTracks, :DocTestSetup, :(using StellarTracks); recursive=true)

# The `format` below makes it so that urls are set to "pretty" if you are pushing them to a hosting service, and basic if you are just using them locally to make browsing easier.

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric) # style=:authoryear
links = InterLinks("BolometricCorrections" => "https://cgarling.github.io/BolometricCorrections.jl/stable/")

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
        "mist.md",
        "api.md",
        "refs.md",
        "doc_index.md"
    ],
    pagesonly=true,
    doctest=false,
    linkcheck=ci, # only check links on CI
    warnonly=[:missing_docs, :linkcheck],
    plugins=[bib, links]
)

deploydocs(;
    repo="github.com/cgarling/StellarTracks.jl",
    versions = ["stable" => "v^", "v#.#"],
    devbranch="main",
    push_preview=true
)

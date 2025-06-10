

function __init__()
    # Register datadep
    register(DataDep("BaSTIv2",
                     """The updated BaSTI stellar models first presented in \
                        Hidalgo et al. 2018.""",
                     "https://github.com/cgarling/StellarTracks.jl/releases/download/v0.0.2/basti_v2.jld2",
                     "b40eca2bb9ad97f0be48966085f2b2f8984ff4eb357e419af1d904ef0c9ac00d"))
end

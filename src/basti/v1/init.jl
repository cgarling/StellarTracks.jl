

function __init__()
    # Register datadep
    register(DataDep("BaSTIv1",
                     """The older BaSTI stellar models presented in \
                        Pietrinferni et al. 2004, 2006, 2013.""",
                     "https://github.com/cgarling/StellarTracks.jl/releases/download/0.0.1/basti_v1.jld2",
                     "13c57c0d11024dedc2ae4634a7ee18ddebae398698dca5359668ac81412cbb84"))
end

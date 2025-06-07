

function __init__()
    # Register datadep
    register(DataDep("BaSTIv1",
                     """The older BaSTI stellar models presented in \
                        Pietrinferni et al. 2004, 2006, 2013.""",
                     "https://github.com/cgarling/StellarTracks.jl/releases/download/0.0.1/basti2013.jld2",
                     "4cf83627031bae97ba63687b1cf172d1e3ed07266b18fbcc81124025245ac264"))
end

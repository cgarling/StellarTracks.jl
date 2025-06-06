

function __init__()
    # Register datadep
    register(DataDep("BaSTIv1",
                     """The older BaSTI stellar models presented in \
                        Pietrinferni et al. 2004, 2006, 2013.""",
                     "https://github.com/cgarling/StellarTracks.jl/releases/download/0.0.1/basti2013.jld2",
                     "7481b0b4e10f673dd47ed0d3046f4c23ac4e91052efadbadf634253ba7c28aef"))
end

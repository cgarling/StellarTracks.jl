using StellarTracks
using Test, SafeTestsets

@testset verbose=true "StellarTracks.jl" begin
    @safetestset "Isochrone Creation" include("isochrones.jl")
end

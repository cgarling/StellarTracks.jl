using StellarTracks
using Test, SafeTestsets

# Run doctests first
using Documenter: DocMeta, doctest
DocMeta.setdocmeta!(StellarTracks, :DocTestSetup, :(using StellarTracks; import BolometricCorrections; using DataDeps: @datadep_str; datadep"PARSECv1.2S"); recursive=true)
doctest(StellarTracks)

@testset verbose=true "StellarTracks.jl" begin
    @safetestset "Isochrone Creation" include("isochrones.jl")
    @safetestset "MIST Tracks" include("mist_tests.jl")
end

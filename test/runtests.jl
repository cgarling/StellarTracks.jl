using StellarTracks
using Test, SafeTestsets

# Run doctests first
using Documenter: DocMeta, doctest
DocMeta.setdocmeta!(StellarTracks, :DocTestSetup, :(using StellarTracks; import BolometricCorrections; using DataDeps: @datadep_str; datadep"PARSECv1.2S"; datadep"MISTv1.2_vvcrit0.0"; datadep"MISTv1.2_vvcrit0.4"; datadep"BaSTIv1"); recursive=true)
doctest(StellarTracks)

@testset verbose=true "StellarTracks.jl" begin
    @safetestset "PARSEC Tracks" include("parsec_tests.jl")
    @safetestset "MIST Tracks" include("mist_tests.jl")
    @safetestset "BaSTIv1 Tracks" include("basti_v1_tests.jl")
end

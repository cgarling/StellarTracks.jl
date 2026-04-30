using StellarTracks
using Test, SafeTestsets

# Run doctests first
using Documenter: DocMeta, doctest
DocMeta.setdocmeta!(StellarTracks, :DocTestSetup, :(using StellarTracks; import BolometricCorrections; using DataDeps: @datadep_str; datadep"PARSECv1.2S"; datadep"MISTv1.2_vvcrit0.0"; datadep"MISTv1.2_vvcrit0.4"; datadep"MISTv2.5_vvcrit0.0_afe_p0.0"; datadep"MISTv2.5_vvcrit0.0_afe_p0.4"; datadep"BaSTIv1"; datadep"BaSTIv2"); recursive=true)
doctest(StellarTracks)

@testset verbose=true "StellarTracks.jl" begin
    @safetestset "PARSEC Tracks" include("parsec_tests.jl")
    @safetestset "MIST Tracks" include("mist_tests.jl")
    @safetestset "BaSTIv1 Tracks" include("basti_v1_tests.jl")
    @safetestset "BaSTIv2 Tracks" include("basti_v2_tests.jl")
    @safetestset "InterpolatedTrack" include("interpolated_track_tests.jl")
    @testset "alpha mass fractions" begin
        @safetestset "BaSTIv1Chemistry (Grevesse & Noels 1993)" include("alpha_mass_fractions/basti_v1_grevesse1993.jl")
        @safetestset "BaSTIv2Chemistry (Caffau+2011)" include("alpha_mass_fractions/basti_v2_caffau2011.jl")
    end
end

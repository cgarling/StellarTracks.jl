using StellarTracks
using StellarTracks: _top_down_indices
using Test

# Retain the high-mass side of AMR folds, comparing ages against the last retained sample.
@testset "Top-down age–mass selection" begin
    @test _top_down_indices([1, 2, 3, 4, 5], [10, 8, 9, 6, 5]) == [5, 4, 3, 1]
    @test _top_down_indices([1, 2, 3, 4, 5], [10, 9, 7, 8, 5]) == [5, 4, 2, 1]
    @test _top_down_indices([1, 2, 3], [10, 8, 8]) == [3, 1]
    @test _top_down_indices([2, 4, 1, 3], [9, 7, 10, 8]) == [2, 4, 1, 3]
    @test _top_down_indices([1, 2, 3], [10, 9, 8]) == [3, 2, 1]
    @test _top_down_indices([1, 2, 3], [10, 10, 10]) == [3]
    @test _top_down_indices(Float64[], Float64[]) == Int[]
    masses = [1.0, 2.0, 3.0, 4.0, 5.0]
    ages = [10.0, 9.0, 7.0, 8.0, 5.0]
    indices = _top_down_indices(masses, ages)
    @test all(<(0), diff(masses[indices]))
    @test all(>(0), diff(ages[indices]))
    @test masses == [1.0, 2.0, 3.0, 4.0, 5.0]
    @test ages == [10.0, 9.0, 7.0, 8.0, 5.0]
end

@testset "Top-down BaSTI mass support" begin
    ts = StellarTracks.BaSTIv1.BaSTIv1TrackSet(0.0003, 0.4, true, false, 0.4)
    @test minimum(mass(ts)) == 0.55f0
    @test ts(minimum(mass(ts))) isa StellarTracks.BaSTIv1.BaSTIv1Track
    @test_throws "Requested mass" ts(0.5)
    @test all(amr -> all(>(0), diff(amr.t)) && all(<(0), diff(amr.u)), ts.AMRs)
end

@testset "_supported_masses drops masses absent from every EEP" begin
    interps = [(t = [1.0, 2.0],), (t = [2.5, 3.0],)]
    @test StellarTracks._supported_masses(interps, [0.5, 1.0, 2.0, 2.5, 3.5]) == [1.0, 2.0, 2.5]
    @test StellarTracks._supported_masses(interps, [1.0, 1.0, 2.0]) == [1.0, 2.0]
end

# Every advertised mass must be usable: top-down selection can drop a mass from every EEP.
@testset "Advertised masses are supported by at least one EEP" begin
    for ts in (StellarTracks.PARSEC.PARSECTrackSet(0.004),
               StellarTracks.MIST.MISTv1TrackSet(0.25, 0.0),
               StellarTracks.MIST.MISTv2TrackSet(-3.5, 0.0, 0.0),
               StellarTracks.BaSTIv1.BaSTIv1TrackSet(0.0003, 0.4, true, false, 0.4),
               StellarTracks.BaSTIv2.BaSTIv2TrackSet(-1.55))
        @test !isempty(mass(ts))
        @test all(m -> any(itp -> first(itp.t) <= m <= last(itp.t), ts.interps), mass(ts))
    end
end

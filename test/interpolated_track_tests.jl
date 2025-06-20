using StellarTracks
using StellarTracks: InterpolatedTrack
using Test
mean(x) = sum(x) / length(x)

for lib in (PARSECLibrary(), MISTLibrary(), BaSTIv1Library(), BaSTIv2Library())
    mean_mh = mean(MH(lib))
    M = 1.51
    track = lib(mean_mh, M)
    @test track isa StellarTracks.InterpolatedTrack
    @test mass(track) ≈ M
    @test MH(track) ≈ mean_mh
    @test Z(track) ≈ Z(chemistry(track), MH(track))
    @test Y(track) ≈ Y(chemistry(track), Z(track))
    @test X(track) ≈ X(chemistry(track), Z(track))
    @test keys(track) == keys(track.track0)
    @test keys(track) == keys(track(9))
    @test extrema(track) isa NTuple{2, <:Number} # NTuple{2, eltype(track)}
    for logage in range(extrema(track)...; length=10)
        itp = track(logage)
        @test itp isa NamedTuple
        @test keys(itp) == keys(track)
    end
end
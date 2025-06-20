using StellarTracks.MIST
using StellarTracks: InterpolatedTrack
using BolometricCorrections.MIST: MISTBCGrid
using TypedTables: Table

using Test

bcg = MISTBCGrid("JWST")

@testset "MIST" begin
    @testset "MISTTrack" begin
        for vvcrit in (0.0, 0.4)
            for feh in MIST.feh_grid
                for M in MIST.mist_massgrid
                    track = MISTTrack(feh, M, vvcrit)
                    @test mass(track) == M
                    @test chemistry(track) == MISTChemistry()
                    extr = extrema(track) # Get logAge limits
                    @test track(extr[1] + (extr[2] - extr[1])/2) isa NamedTuple
                    @test MH(track) == feh
                    @test Z(track) == Z(chemistry(track), MH(track))
                    @test Y(track) == Y(chemistry(track), Z(track))
                    @test X(track) == X(chemistry(track), Z(track))
                    @test post_rgb(track) isa Bool
                    @test eltype(track) == eltype(MIST.feh_grid)
                end
            end
        end
    end
    @testset "MISTTrackSet" begin
        for vvcrit in (0.0, 0.4)
            for feh in MIST.feh_grid
                trackset = MISTTrackSet(feh, vvcrit)
                @test trackset(1.0) isa MISTTrack
                @test mass(trackset) == MIST.mist_massgrid
                @test chemistry(trackset) == MISTChemistry()
                @test MH(trackset) == feh
                @test Z(trackset) == Z(chemistry(trackset), MH(trackset))
                @test Y(trackset) == Y(chemistry(trackset), Z(trackset))
                @test X(trackset) == X(chemistry(trackset), Z(trackset))
                @test post_rgb(trackset) == true
                @test eltype(trackset) == eltype(MIST.feh_grid)

                # Test isochrone, with and without BCs
                @test isochrone(trackset, 10.0) isa NamedTuple
                # Interpolate MISTBCGrid to appropriate feh, no reddening
                bct = bcg(feh, 0.0)
                @test isochrone(trackset, bct, 10.0) isa Table
            end
        end
    end
    @testset "MISTLibrary" begin
        for vvcrit in (0.0, 0.4)
            tracklib = MISTLibrary(vvcrit)
            @test_throws "Not yet implemented." tracklib(-4)
            @test tracklib(-2.05, 1.05) isa InterpolatedTrack
            @test chemistry(tracklib) == MISTChemistry()
            @test MH(tracklib) == MIST.feh_grid
            @test Z(tracklib) == Z.(chemistry(tracklib), MH(tracklib))
            @test Y(tracklib) == Y.(chemistry(tracklib), Z(tracklib))
            @test X(tracklib) == X.(chemistry(tracklib), Z(tracklib))
            @test post_rgb(tracklib) == true
            @test eltype(tracklib) == eltype(MIST.feh_grid)

            # Test isochrone, with and without BCs
            @test isochrone(tracklib, 10.0, -2.15) isa NamedTuple
            # Test passing full stellar track library with single BC table
            iso1 = isochrone(tracklib, bcg(-2.15, 0.0), 10.0, -2.15)
            @test iso1 isa Table
            iso2 = isochrone(tracklib, bcg, 10.0, -2.15, 0.0)
            @test iso2 isa Table
            # These should both be equivalent
            @test iso1 == iso2
            # Test that we can interpolate across the range of valid MH
            for feh in range(extrema(MH(tracklib))...; length=10)
                @test_nowarn isochrone(tracklib, 10.0, feh)
                @test_nowarn isochrone(tracklib, bcg, 10.0, feh, 0.0)
            end
            # Make sure error is thrown for out of MH range
            @test_throws DomainError isochrone(tracklib, 10.0, -5)
            @test_throws DomainError isochrone(tracklib, 10.0, 1.0)
        end
    end
    
end

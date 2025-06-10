using StellarTracks.PARSEC
using BolometricCorrections.MIST: MISTBCGrid
using TypedTables: Table

using Test

bcg = MISTBCGrid("JWST")

@testset "PARSEC" begin
    @testset "PARSECTrack + PARSECTrackSet" begin
        # Test PARSECTrack() constructor; inefficient as currently implemented
        # because it requires loading an entire PARSECTrackSet, but convenient
        @test PARSECTrack(0.0001, 1.05) isa PARSECTrack
        for z in PARSEC.zgrid
            trackset = PARSECTrackSet(z)
            @test trackset(1.0) isa PARSECTrack
            # @test mass(trackset) == PARSEC.parsec_massgrid
            @test chemistry(trackset) == PARSECChemistry()
            @test Z(trackset) == z
            @test MH(trackset) == MH(chemistry(trackset), Z(trackset))
            @test Y(trackset) == Y(chemistry(trackset), Z(trackset))
            @test X(trackset) == X(chemistry(trackset), Z(trackset))
            @test post_rgb(trackset) == true
            @test eltype(trackset) == eltype(PARSEC.track_type)

            # Test isochrone, with and without BCs
            @test isochrone(trackset, 10.0) isa NamedTuple
            # Interpolate MISTBCGrid to appropriate feh, no reddening
            bct = bcg(MH(trackset), 0.0)
            @test isochrone(trackset, bct, 10.0) isa Table
            for M in range(extrema(mass(trackset))...; length=10)
                track = trackset(M)
                @test mass(track) == M
                @test chemistry(track) == PARSECChemistry()
                extr = extrema(track) # Get logAge limits
                @test track(extr[1] + (extr[2] - extr[1])/2) isa NamedTuple
                @test Z(track) == z
                @test MH(track) == MH(chemistry(track), Z(track))
                @test Y(track) == Y(chemistry(track), Z(track))
                @test X(track) == X(chemistry(track), Z(track))
                @test post_rgb(track) isa Bool
                @test eltype(track) == PARSEC.track_type
            end
        end
    end

    @testset "PARSECLibrary" begin
        tracklib = PARSECLibrary()
        @test_throws "Not yet implemented." tracklib(0.0001234)
        @test chemistry(tracklib) == PARSECChemistry()
        @test Z(tracklib) == PARSEC.zgrid
        @test MH(tracklib) == MH.(chemistry(tracklib), Z(tracklib))
        @test Y(tracklib) == Y.(chemistry(tracklib), Z(tracklib))
        @test X(tracklib) == X.(chemistry(tracklib), Z(tracklib))
        @test post_rgb(tracklib) == true
        @test eltype(tracklib) == PARSEC.track_type

        # Test isochrone, with and without BCs
        @test isochrone(tracklib, 10.0, -1.234) isa NamedTuple
        # Test passing full stellar track library with single BC table
        iso1 = isochrone(tracklib,
                         bcg(MH(chemistry(bcg), Z(chemistry(tracklib), -1.234)), 0.0),
                         10.0, -1.234)
        @test iso1 isa Table
        iso2 = isochrone(tracklib, bcg, 10.0, -1.234, 0.0)
        @test iso2 isa Table
        # These should both be equivalent
        @test iso1 == iso2
        # Test that we can interpolate across the range of valid MH
        for z in range(extrema(Z(tracklib))...; length=10)
            @test_nowarn isochrone(tracklib, 10.0, z)
            @test_nowarn isochrone(tracklib, bcg, 10.0, z, 0.0)
        end
        # Make sure error is thrown for out of MH range
        @test_throws DomainError isochrone(tracklib, 10.0, -5)
        @test_throws DomainError isochrone(tracklib, 10.0, 1.0)
    end
    
end

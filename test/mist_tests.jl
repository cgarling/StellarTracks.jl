using StellarTracks.MIST
using StellarTracks: InterpolatedTrack
using BolometricCorrections.MIST: MISTv1BCGrid, MISTv2BCGrid
using TypedTables: Table

using Test

@test gridname(MISTv1Track) isa String
@test gridname(MISTv1TrackSet) isa String
@test gridname(MISTv1Library) isa String

bcgv1 = MISTv1BCGrid("JWST")
bcgv2 = MISTv2BCGrid("JWST")

@testset "MIST" begin
    @testset "MISTv1Track" begin
        for vvcrit in (0.0, 0.4)
            for feh in MIST.feh_grid
                for M in MIST.mass_grid
                    track = MISTv1Track(feh, M, vvcrit)
                    @test gridname(track) isa String
                    @test mass(track) == M
                    @test chemistry(track) == MISTv1Chemistry()
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
    @testset "MISTv1TrackSet" begin
        for vvcrit in (0.0, 0.4)
            for feh in MIST.feh_grid
                trackset = MISTv1TrackSet(feh, vvcrit)
                @test gridname(trackset) isa String
                @test trackset(1.0) isa MISTv1Track
                @test mass(trackset) == MIST.mass_grid
                @test chemistry(trackset) == MISTv1Chemistry()
                @test MH(trackset) == feh
                @test Z(trackset) == Z(chemistry(trackset), MH(trackset))
                @test Y(trackset) == Y(chemistry(trackset), Z(trackset))
                @test X(trackset) == X(chemistry(trackset), Z(trackset))
                @test post_rgb(trackset) == true
                @test eltype(trackset) == eltype(MIST.feh_grid)

                # Test isochrone, with and without BCs
                @test isochrone(trackset, 10.0) isa NamedTuple
                # Interpolate MISTv1BCGrid to appropriate feh, no reddening
                bctv1 = bcgv1(feh, 0.0)
                @test isochrone(trackset, bctv1, 10.0) isa Table
                if extrema(bcgv2).feh[1] <= feh <= extrema(bcgv2).feh[2]
                    # Interpolate MISTv2BCGrid to appropriate feh, no reddening, solar α-abundance
                    bctv2 = bcgv2(feh, 0.0, 0.0)
                    @test isochrone(trackset, bctv2, 10.0) isa Table
                end
            end
        end
    end
    @testset "MISTv1Library" begin
        for vvcrit in (0.0, 0.4)
            tracklib = MISTv1Library(vvcrit)
            @test gridname(tracklib) isa String
            @test tracklib(-2.05, 1.05) isa InterpolatedTrack
            @test chemistry(tracklib) == MISTv1Chemistry()
            @test MH(tracklib) == MIST.feh_grid
            @test Z(tracklib) == Z.(chemistry(tracklib), MH(tracklib))
            @test Y(tracklib) == Y.(chemistry(tracklib), Z(tracklib))
            @test X(tracklib) == X.(chemistry(tracklib), Z(tracklib))
            @test post_rgb(tracklib) == true
            @test eltype(tracklib) == eltype(MIST.feh_grid)

            # Test isochrone, with and without BCs
            @test isochrone(tracklib, 10.0, -2.15) isa NamedTuple
            # Test passing full stellar track library with single BC table
            iso1 = isochrone(tracklib, bcgv1(-2.15, 0.0), 10.0, -2.15)
            @test iso1 isa Table
            iso2 = isochrone(tracklib, bcgv1, 10.0, -2.15, 0.0)
            @test iso2 isa Table
            # These should both be equivalent
            @test iso1 == iso2

            iso1_v2 = isochrone(tracklib, bcgv2(-2.15, 0.0, 0.0), 10.0, -2.15)
            @test iso1_v2 isa Table
            iso2_v2 = isochrone(tracklib, bcgv2, 10.0, -2.15, 0.0, 0.0)
            @test iso2_v2 isa Table
            # These will not be the same because the iso_v2 call will convert the tracklib's feh to Z and then back to feh for the BC interpolation,
            # which will not be exactly the same as the original feh due to differences in the solar abundance pattern between the tracklib and BC grid.
            # @test iso1_v2 == iso2_v2

            # Test that we can interpolate across the range of valid MH
            for feh in range(extrema(MH(tracklib))...; length=10)
                @test_nowarn isochrone(tracklib, 10.0, feh)
                @test_nowarn isochrone(tracklib, bcgv1, 10.0, feh, 0.0)
                feh_v2 = MH(chemistry(bcgv2), Z(chemistry(tracklib), feh))
                if extrema(bcgv2).feh[1] <= feh_v2 <= extrema(bcgv2).feh[2]
                    @test_nowarn isochrone(tracklib, bcgv2, 10.0, feh, 0.0, 0.0)
                end
            end
            # Make sure error is thrown for out of MH range
            @test_throws DomainError isochrone(tracklib, 10.0, -5)
            @test_throws DomainError isochrone(tracklib, 10.0, 1.0)
        end
    end
    
end

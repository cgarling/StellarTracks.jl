using StellarTracks.BaSTIv1
using StellarTracks: InterpolatedTrack
using BolometricCorrections.MIST: MISTBCGrid
using TypedTables: Table

using Test

bcg = MISTBCGrid("JWST")

@testset "BaSTIv1" begin
    @testset "BaSTIv1Track + BaSTIv1TrackSet" begin
        # Test BaSTIv1Track() constructor; inefficient as currently implemented
        # because it requires loading an entire BaSTIv1TrackSet, but convenient
        @test BaSTIv1Track(0.0001, 1.2, 0.0, true, false, 0.4) isa BaSTIv1Track
        for agb in (true, false)
            for η in BaSTIv1.ηgrid
                # First and last values of BaSTIv1.zgrid only available
                # for η = 0.4, without AGB extension
                if ~agb && η ≈ 0.4
                    zg = BaSTIv1.zgrid
                else
                    zg = BaSTIv1.zgrid[begin+1:end-1]
                end
                for canonical in (true, false)
                    for α_fe in BaSTIv1.αFegrid
                        for z in zg
                            # println(canonical, " ", α_fe, " ", z)
                            trackset = BaSTIv1TrackSet(z, α_fe, canonical, agb, η)
                            @test trackset(minimum(mass(trackset))+0.1) isa BaSTIv1Track
                            # @test mass(trackset) == PARSEC.parsec_massgrid
                            @test chemistry(trackset) == BaSTIv1Chemistry()
                            @test Z(trackset) == z
                            @test MH(trackset) == MH(chemistry(trackset), Z(trackset))
                            @test Y(trackset) == Y(chemistry(trackset), Z(trackset))
                            @test X(trackset) == X(chemistry(trackset), Z(trackset))
                            @test post_rgb(trackset) == true
                            @test eltype(trackset) == eltype(BaSTIv1.track_type)

                            # Test isochrone, with and without BCs
                            @test isochrone(trackset, 9.0) isa NamedTuple
                            # Interpolate MISTBCGrid to appropriate feh, no reddening
                            bct = bcg(MH(trackset), 0.0)
                            @test isochrone(trackset, bct, 9.0) isa Table
                            for M in range(extrema(mass(trackset))...; length=10)
                                # println(M)
                                track = trackset(M)
                                @test mass(track) == M
                                @test chemistry(track) == BaSTIv1Chemistry()
                                extr = extrema(track) # Get logAge limits
                                @test track(extr[1] + (extr[2] - extr[1])/2) isa NamedTuple
                                @test Z(track) == z
                                @test MH(track) == MH(chemistry(track), Z(track))
                                @test Y(track) == Y(chemistry(track), Z(track))
                                @test X(track) == X(chemistry(track), Z(track))
                                @test post_rgb(track) isa Bool
                                @test eltype(track) == BaSTIv1.track_type
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "BaSTIv1Library" begin
        for agb in (true, false)
            for η in BaSTIv1.ηgrid
                # First and last values of BaSTIv1.zgrid only available
                # for η = 0.4, without AGB extension
                if ~agb && η ≈ 0.4
                    zg = BaSTIv1.zgrid
                else
                    zg = BaSTIv1.zgrid[begin+1:end-1]
                end
                for canonical in (true, false)
                    for α_fe in BaSTIv1.αFegrid
                        tracklib = BaSTIv1Library(α_fe, canonical, agb, η)
                        tracklib(-1.05, 1.51) isa InterpolatedTrack
                        @test chemistry(tracklib) == BaSTIv1Chemistry()
                        @test Z(tracklib) == zg # BaSTIv1.zgrid
                        @test MH(tracklib) == MH.(chemistry(tracklib), Z(tracklib))
                        @test Y(tracklib) == Y.(chemistry(tracklib), Z(tracklib))
                        @test X(tracklib) == X.(chemistry(tracklib), Z(tracklib))
                        @test post_rgb(tracklib) == true
                        @test eltype(tracklib) == BaSTIv1.track_type

                        # Test isochrone, with and without BCs
                        @test isochrone(tracklib, 9.0, -1.234) isa NamedTuple
                        # Test passing full stellar track library with single BC table
                        iso1 = isochrone(tracklib,
                                         bcg(MH(chemistry(bcg), Z(chemistry(tracklib), -1.234)), 0.0),
                                         9.0, -1.234)
                        @test iso1 isa Table
                        iso2 = isochrone(tracklib, bcg, 9.0, -1.234, 0.0)
                        @test iso2 isa Table
                        # These should both be equivalent
                        @test iso1 == iso2
                        # Test that we can interpolate across the range of valid MH
                        for mh in range(extrema(MH(tracklib))...; length=10)
                            @test_nowarn isochrone(tracklib, 9.0, mh)
                            @test_nowarn isochrone(tracklib, bcg, 9.0, mh, 0.0)
                        end
                        # Make sure error is thrown for out of MH range
                        @test_throws DomainError isochrone(tracklib, 9.0, minimum(MH(tracklib)) - 1)
                        @test_throws DomainError isochrone(tracklib, 9.0, maximum(MH(tracklib)) + 1)
                    end
                end
            end
        end
    end
end

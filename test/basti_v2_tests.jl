using StellarTracks.BaSTIv2
using StellarTracks.BaSTIv2: track_type # number type that the data is represented as
using BolometricCorrections.MIST: MISTBCGrid
using TypedTables: Table

using Test

bcg = MISTBCGrid("JWST")

param_sets = ((α=0.0, canonical = false, diffusion = true, yp = 0.247, η = 0.3),
              (α=0.0, canonical = false, diffusion = false, yp = 0.247, η = 0.0),
              (α=0.0, canonical = true, diffusion = false, yp = 0.247, η = 0.0),
              (α=0.0, canonical = false, diffusion = false, yp = 0.247, η = 0.3),
              (α=-0.2, canonical = false, diffusion = true, yp = 0.247, η = 0.3))

@testset "BaSTIv2" begin
    @testset "BaSTIv2Track + BaSTIv2TrackSet" begin
        # Test BaSTIv2Track() constructor; inefficient as currently implemented
        # because it requires loading an entire BaSTIv1TrackSet, but convenient
        @test BaSTIv2Track(-2.2, 1.05, 0.0, false, true, 0.247, 0.3) isa BaSTIv2Track
        for p in param_sets
            vp = values(p)
            for feh in BaSTIv2.feh_grid
                chem = BaSTIv2Chemistry(convert(track_type, p.α),
                                        convert(track_type, p.yp))
                # println(canonical, " ", α_fe, " ", z)
                trackset = BaSTIv2TrackSet(feh, vp...)
                @test trackset(minimum(mass(trackset))+0.1) isa BaSTIv2Track
                # @test mass(trackset) == PARSEC.parsec_massgrid
                @test chemistry(trackset) == chem
                # @test Z(trackset) == z
                @test isapprox(MH(trackset), MH(chemistry(trackset), Z(trackset)); atol=1e-5)
                @test Y(trackset) ≈ Y(chemistry(trackset), Z(trackset))
                @test X(trackset) ≈ X(chemistry(trackset), Z(trackset))
                @test post_rgb(trackset) == true
                @test eltype(trackset) == eltype(BaSTIv2.track_type)

                # Test isochrone, with and without BCs
                @test isochrone(trackset, 9.0) isa NamedTuple
                # Interpolate MISTBCGrid to appropriate feh, no reddening
                bct = bcg(MH(trackset), 0.0)
                @test isochrone(trackset, bct, 9.0) isa Table
                for M in range(extrema(mass(trackset))...; length=10)
                    # println(M)
                    track = trackset(M)
                    @test mass(track) == M
                    @test chemistry(track) == chem
                    extr = extrema(track) # Get logAge limits
                    @test track(extr[1] + (extr[2] - extr[1])/2) isa NamedTuple
                    # @test Z(track) == z
                    @test isapprox(MH(track), MH(chemistry(track), Z(track)); atol=1e-5)
                    @test Y(track) ≈ Y(chemistry(track), Z(track))
                    @test X(track) ≈ X(chemistry(track), Z(track))
                    @test post_rgb(track) isa Bool
                    @test eltype(track) == BaSTIv2.track_type
                end
            end
        end
    end

    @testset "BaSTIv2Library" begin
        for p in param_sets
            vp = values(p)
            chem = BaSTIv2Chemistry(convert(track_type, p.α),
                                    convert(track_type, p.yp))
            
            tracklib = BaSTIv2Library(vp...)
            @test_throws "Not yet implemented." tracklib(0.0001234)
            @test chemistry(tracklib) == chem
            # @test Z(tracklib) == BaSTIv1.zgrid
            @test MH(tracklib) ≈ MH.(chemistry(tracklib), Z(tracklib))
            @test Y(tracklib) ≈ Y.(chemistry(tracklib), Z(tracklib))
            @test X(tracklib) ≈ X.(chemistry(tracklib), Z(tracklib))
            @test post_rgb(tracklib) == true
            @test eltype(tracklib) == BaSTIv2.track_type

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
            for z in range(extrema(Z(tracklib))...; length=10)
                @test_nowarn isochrone(tracklib, 9.0, z)
                @test_nowarn isochrone(tracklib, bcg, 9.0, z, 0.0)
            end
            # Make sure error is thrown for out of MH range
            @test_throws DomainError isochrone(tracklib, 9.0, -5)
            @test_throws DomainError isochrone(tracklib, 9.0, 1.0)
        end
    end
end

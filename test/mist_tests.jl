using StellarTracks.MIST
using StellarTracks: InterpolatedTrack
using BolometricCorrections.MIST: MISTv1BCGrid, MISTv2BCGrid
using TypedTables: Table

using Test

@test gridname(MISTv1Track) isa String
@test gridname(MISTv1TrackSet) isa String
@test gridname(MISTv1Library) isa String
@test gridname(MISTv2Track) isa String
@test gridname(MISTv2TrackSet) isa String
@test gridname(MISTv2Library) isa String

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
            iso2_v2 = isochrone(tracklib, bcgv2, 10.0, -2.15, 0.0)
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
                    @test_nowarn isochrone(tracklib, bcgv2, 10.0, feh, 0.0)
                end
            end
            # Make sure error is thrown for out of MH range
            @test_throws DomainError isochrone(tracklib, 10.0, -5)
            @test_throws DomainError isochrone(tracklib, 10.0, 1.0)
        end
    end

    @testset "MISTv2Track" begin
        # Spot-check a representative subset of (vvcrit, afe, feh, mass) combinations
        for (vvcrit, afe) in ((0.0, 0.0), (0.4, 0.0), (0.4, 0.4))
            for feh in [-2.0, 0.0]  # subset of feh_grid_v2
                for M in [0.1, 0.81, 2.0]  # subset of mass_grid_v2
                    track = MISTv2Track(feh, M, vvcrit, afe)
                    @test gridname(track) isa String
                    @test mass(track) == M
                    @test chemistry(track) == MISTv2Chemistry()
                    @test alphaFe(track) == afe
                    @test FeH(track) == feh
                    extr = extrema(track)
                    @test track(extr[1] + (extr[2] - extr[1])/2) isa NamedTuple
                    @test isapprox(Z(track), Z(chemistry(track), MH(track)); rtol=1e-5)
                    @test isapprox(Y(track), Y(chemistry(track), Z(track)); rtol=1e-5)
                    @test isapprox(X(track), X(chemistry(track), Z(track)); rtol=1e-5)
                    @test post_rgb(track) isa Bool
                    @test eltype(track) == eltype(MIST.feh_grid_v2)
                end
            end
        end
    end

    @testset "MISTv2TrackSet" begin
        for (vvcrit, afe) in ((0.0, 0.0), (0.0, 0.4), (0.4, 0.4))
            for feh in MIST.feh_grid_v2_for(afe)
                trackset = MISTv2TrackSet(feh, vvcrit, afe)
                @test gridname(trackset) isa String
                @test trackset(minimum(mass(trackset)) + 0.1) isa MISTv2Track
                @test chemistry(trackset) == MISTv2Chemistry()
                @test alphaFe(trackset) == afe
                @test FeH(trackset) == feh
                @test isapprox(MH(trackset), MH(chemistry(trackset), Z(trackset)); atol=1e-5)
                @test isapprox(Z(trackset), Z(chemistry(trackset), MH(trackset)); atol=1e-5)
                @test isapprox(Y(trackset), Y(chemistry(trackset), Z(trackset)); rtol=1e-5)
                @test isapprox(X(trackset), X(chemistry(trackset), Z(trackset)); rtol=1e-5)
                @test post_rgb(trackset) == true
                @test eltype(trackset) == eltype(MIST.feh_grid_v2)
                @test isochrone(trackset, 10.0) isa NamedTuple
                # Test isochrone with MISTv2BCGrid
                if extrema(bcgv2).feh[1] <= FeH(trackset) <= extrema(bcgv2).feh[2]
                    bctv2 = bcgv2(FeH(trackset), afe, 0.0)
                    @test isochrone(trackset, bctv2, 10.0) isa Table
                end
            end
        end
    end

    @testset "MISTv2Library" begin
        for (vvcrit, afe) in ((0.0, 0.0), (0.0, 0.4), (0.4, 0.4))
            tracklib = MISTv2Library(vvcrit, afe)
            @test gridname(tracklib) isa String
            @test tracklib(-2.05, 1.05) isa InterpolatedTrack
            @test chemistry(tracklib) == MISTv2Chemistry()
            @test alphaFe(tracklib) == afe
            @test FeH(tracklib) == MIST.feh_grid_v2_for(afe)
            @test MH(tracklib) ≈ MH.(chemistry(tracklib), Z(tracklib))
            @test Y(tracklib) ≈ Y.(chemistry(tracklib), Z(tracklib))
            @test X(tracklib) ≈ X.(chemistry(tracklib), Z(tracklib))
            @test post_rgb(tracklib) == true
            @test eltype(tracklib) == eltype(MIST.feh_grid_v2)

            @test isochrone(tracklib, 10.0, -2.15) isa NamedTuple

            # Test BC integration with MISTv2BCGrid.
            # MISTv2Library and MISTv2BCGrid both use MISTv2Chemistry, so bc_mh ≈ mh
            # and afe is inferred from alphaFe(tracklib). The two call paths should
            # produce identical results.
            mh_test = eltype(tracklib)(-2.15)
            bc_mh = MH(chemistry(bcgv2), Z(chemistry(tracklib), mh_test))
            # MISTv2BCGrid (bcgv2) expects feh and afe, so convert from mh and tracklib's
            # chemistry to get the correct feh for the BC interpolation
            bc_feh = FeH(chemistry(bcgv2), bc_mh, alphaFe(tracklib))
            iso1_v2 = isochrone(tracklib, bcgv2(bc_feh, afe, 0.0), 10.0, mh_test)
            @test iso1_v2 isa Table
            iso2_v2 = isochrone(tracklib, bcgv2, 10.0, mh_test, 0.0)
            @test iso2_v2 isa Table
            @test iso1_v2 == iso2_v2

            # Test that we can interpolate across the range of valid MH
            for mh in range(extrema(MH(tracklib))...; length=10)
                @test_nowarn isochrone(tracklib, 10.0, mh)
                # bc_mh == mh for MISTv2 (same chemistry), but check it's within
                # the BC grid's valid [Fe/H] range before calling
                bc_mh = MH(chemistry(bcgv2), Z(chemistry(tracklib), mh))
                bc_feh = FeH(chemistry(bcgv2), bc_mh, alphaFe(tracklib))
                if extrema(bcgv2).feh[1] <= bc_feh <= extrema(bcgv2).feh[2]
                    @test_nowarn isochrone(tracklib, bcgv2, 10.0, mh, 0.0)
                end
            end
            # Make sure error is thrown for out of MH range
            @test_throws DomainError isochrone(tracklib, 10.0, -5)
            @test_throws DomainError isochrone(tracklib, 10.0, 1.0)
        end
    end

end

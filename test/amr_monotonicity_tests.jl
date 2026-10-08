using StellarTracks
using StellarTracks: _amr, _increasing_subsequence
using StellarTracks.PARSEC: PARSECTrackSet, PARSECLibrary, eep_idxs, zgrid, ygrid, excluded_masses, _fix_post_rgb_eeps
using StellarTracks.MIST: MISTv1TrackSet
using DataDeps: @datadep_str
import JLD2
using Test

# Shared knots: if every track is at least as old at EEP e+1 as at EEP e, the e+1 mass is never lower at any age
@testset "Linear log-mass AMRs preserve EEP order" begin
    m = [5.0, 3.0, 2.0, 1.2, 0.8]
    t0 = [7.6, 8.4, 8.9, 9.5, 10.0]
    t1 = t0 .+ [0.001, 0.2, 0.0, 0.05, 0.3]
    a0, a1 = _amr(m, t0), _amr(m, t1)
    ages = range(maximum(first.((t0, t1))), minimum(last.((t0, t1))); length=1001)
    @test all(exp10(a1(la)) >= exp10(a0(la)) for la in ages)
    @test exp10(a0(t0[3])) ≈ m[3]
end

@testset "Longest increasing subsequence" begin
    @test _increasing_subsequence([1, 3, 2, 4, 5]) == [1, 3, 4, 5]
    # A small decrease before a long increasing run keeps the run rather than the point before the decrease
    x = [1.0, 2.0, 3.0, 2.5, 2.6, 2.7, 2.8, 2.9]
    @test x[_increasing_subsequence(x)] == [1.0, 2.0, 2.5, 2.6, 2.7, 2.8, 2.9]
    @test length(_increasing_subsequence([2.0, 2.0, 2.0])) == 1
    @test isempty(_increasing_subsequence(Float64[]))
end

@testset "Isochrone initial masses increase" begin
    p = PARSECLibrary()
    for la in 6.8:0.1:10.1, mh in (-2.0, -1.19, -0.75, -0.02, 0.3)
        iso = isochrone(p, la, mh)
        @test all(>(0), diff(iso.m_ini))
    end
    ts = MISTv1TrackSet(0.25, 0.0)
    for la in 6.8:0.2:10.1
        @test all(>(0), diff(isochrone(ts, la).m_ini))
    end
    # At [Fe/H] = -3 most tracks near 2 Msun end before the TP-AGB, so the isochrone mass dips slightly at its start
    # (EEP 809); keeping only points above the running maximum removed all 601 TP-AGB points at logAge 8.8
    @test count(e -> 809 <= e <= 1409, isochrone(MISTv1TrackSet(-3.0, 0.0), 8.8).eep) > 500
end

# A track interpolated at an isochrone point's mass passes through that point: the track's age at the EEP is the exact
# inverse of the AMR that gave the isochrone mass
@testset "Track interpolation inverts the AMR" begin
    ts = PARSECTrackSet(0.004)
    iso = isochrone(ts, 9.0)
    for j in (100, 400, 700)
        track = ts(iso.m_ini[j])
        k = argmin(abs.(track.data.logAge .- 9.0))
        @test track.data.logAge[k] ≈ 9.0 atol=1e-8
        @test track.data.logTe[k] ≈ iso.logTe[j]
    end
end

@testset "PARSEC excluded tracks" begin
    for (z, y, m) in zip(zgrid, ygrid, excluded_masses)
        raw = JLD2.load_object(joinpath(datadep"PARSECv1.2S", "Z$(z)Y$(y).jld2"))
        @test m in raw.m_ini # fails if the list no longer matches the cached grid masses
    end
    ts = PARSECTrackSet(0.02)
    @test !(0.25 in mass(ts))
    # With the 0.25 Msun track, the isochrone mass at logAge 9.45 dropped from 0.218 to 0.138 between EEPs 265 and 266
    m265, m266 = (exp10(ts.AMRs[findfirst(==(e), ts.eeps)](9.45)) for e in (265, 266))
    @test isapprox(m266, m265; rtol=0.05)
end

@testset "PARSEC post-RGB EEP corrections" begin
    raw = JLD2.load_object(joinpath(datadep"PARSECv1.2S", "Z0.014Y0.273.jld2"))
    raw_eep = copy(raw.eep)
    fixed = _fix_post_rgb_eeps(raw)
    @test raw.eep == raw_eep # input not mutated
    track(tab, m) = (t = tab[tab.m_ini .== m]; t[sortperm(t.eep)])
    # 0.76 Msun has only HB-file rows; 1.9 Msun joins MS and HB files; 2.3 Msun is a full track
    @test minimum(track(fixed, 0.76).eep) == eep_idxs.HE_BEG
    hb = track(fixed, 1.9)
    @test !any(e -> eep_idxs.RG_TIP <= e < eep_idxs.HE_BEG, hb.eep)
    @test maximum(hb.eep) == eep_idxs.TPAGB_BEG - 1
    @test track(fixed, 2.3).eep == track(raw, 2.3).eep
    # Fraction of the core-He-burning lifetime reached at a mid-phase EEP; before the correction it is ~0.01 for the
    # HB-file track (1.9 Msun) and ~0.94 for the full track (2.0 Msun)
    function frac(t, e)
        a = exp10.(t.logAge)
        a0, a1 = a[t.eep .== eep_idxs.HE_BEG][1], a[t.eep .== eep_idxs.END_CHEB][1]
        return (a[t.eep .== e][1] - a0) / (a1 - a0)
    end
    @test abs(frac(track(fixed, 1.9), 1381) - frac(track(fixed, 2.0), 1381)) < 0.05
    for m in (0.76, 1.9, 2.3)
        @test issorted(track(fixed, m).logAge)
    end
    # Phase end points keep their ages and properties
    f, r = track(fixed, 2.3), track(raw, 2.3)
    for e in (eep_idxs.HE_BEG, eep_idxs.END_CHEB, eep_idxs.TPAGB_BEG - 1)
        @test f.logAge[f.eep .== e] ≈ r.logAge[r.eep .== e]
        @test f.Mbol[f.eep .== e] ≈ r.Mbol[r.eep .== e]
    end
end

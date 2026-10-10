# Precompile workload: small synthetic track sets with the column types of each library's data files, written to and
# read back from temporary JLD2 files the way the constructors read them, so that first calls to the constructors,
# `isochrone`, and track interpolation do not compile. The data files themselves are downloaded on first use, so they
# cannot be used here.
using PrecompileTools: @compile_workload
import JLD2

# Runs the workload. Downstream packages whose dependencies invalidate this package's cached code (e.g.,
# LoopVectorization.jl, through StarFormationHistories.jl) can call it in their own precompile workloads, so that their
# caches hold this code compiled with those dependencies loaded.
function _precompile_workload()
    masses = [0.8, 1.0, 1.5, 2.0, 3.0]
    # Ages increase along each track and decrease with mass at each EEP, so every EEP has an AMR over all tracks
    age(m, e, n) = 1e10 * m^-2.5 * e / n
    rows(n) = [(m, e) for m in masses for e in 1:n]
    mist_raw(m, n=30) = Table(star_age = [age(m, e, n) for e in 1:n], log_L = [log10(m) + 0.01e for e in 1:n],
                              log_Teff = fill(3.7, n), log_g = fill(4.4, n), log_surf_cell_z = fill(-2.0, n))
    # Every EEP through TPAGB_BEG, so that `_fix_post_rgb_eeps` resamples both post-RGB segments
    function parsec_raw(n=PARSEC.eep_idxs.TPAGB_BEG - 1)
        r = rows(n)
        return Table(logAge = [log10(age(m, e, n)) for (m, e) in r], mass = first.(r), logTe = [3.7 + 1e-4e for (m, e) in r],
                     Mbol = [4.7 - 2.5log10(m) - 1e-3e for (m, e) in r], logg = fill(4.4, length(r)), C_O = fill(0.3, length(r)),
                     m_ini = first.(r), eep = last.(r))
    end
    function basti_raw(v, n=30)
        r = rows(n)
        eep, m_ini = Int16.(last.(r)), Float32.(first.(r))
        logAge, logL, logTe = [log10(age(m, e, n)) for (m, e) in r], Float32[log10(m) + 0.01e for (m, e) in r], fill(3.7f0, length(r))
        return v == 1 ? Table(; eep, m_ini, mass = m_ini, logAge, logL, logTe) : Table(; eep, m_ini, logAge, logL, logTe, mass = m_ini)
    end
    function exercise(lib)
        ts = first(lib.ts)
        isochrone(ts, 9.5)
        ts(1.2)
        mh = (MH(lib)[1] + MH(lib)[2]) / 2 # between the first two track sets
        isochrone(lib, 9.5, mh)
        lib(mh, 1.2)
    end
    mktempdir() do dir
        mistfiles = [joinpath(dir, "mist$i.jld2") for i in eachindex(masses)]
        foreach((f, m) -> JLD2.save_object(f, mist_raw(m)), mistfiles, masses)
        parsecfile = joinpath(dir, "parsec.jld2")
        JLD2.save_object(parsecfile, parsec_raw())
        bastifiles = [joinpath(dir, "basti$v.jld2") for v in 1:2]
        for v in 1:2
            JLD2.jldopen(f -> (f["a/b"] = basti_raw(v)), bastifiles[v], "w")
        end
        data = MIST._mist_table(mistfiles, masses)
        ts = [MIST.MISTv1TrackSet(data, feh, 0.0) for feh in (-1.0, -0.5)]
        # Library types follow the argument types (vvcrit, α_fe), so cover both Int and Float64
        foreach(vvcrit -> exercise(MIST.MISTv1Library(ts, [-1.0, -0.5], vvcrit)), (0, 0.0))
        exercise(MIST.MISTv2Library([MIST.MISTv2TrackSet(data, feh, 0.0, 0.0) for feh in (-1.0, -0.5)], [-1.0, -0.5], 0.0, 0.0))
        # PARSECLibrary takes its metallicities from the full zgrid, so the track set vector has that length
        t = PARSEC._parsec_table(parsecfile, 0.0)
        ts = [PARSEC.PARSECTrackSet(t, PARSEC.zgrid[1]); fill(PARSEC.PARSECTrackSet(t, PARSEC.zgrid[2]), length(PARSEC.zgrid) - 1)]
        exercise(PARSEC.PARSECLibrary(ts))
        d = JLD2.load(bastifiles[1], "a/b")
        ts = [BaSTIv1.BaSTIv1TrackSet(d, z, 0.0, true, false, 0.4) for z in (0.001, 0.002)]
        foreach(α_fe -> exercise(BaSTIv1.BaSTIv1Library(ts, (Z = [0.001, 0.002], α_fe = α_fe, canonical = true, agb = false, η = 0.4))), (0, 0.0))
        d = JLD2.load(bastifiles[2], "a/b")
        ts = [BaSTIv2.BaSTIv2TrackSet(d, feh, 0.0, false, true, 0.247, 0.3) for feh in (-1.05, -0.9)]
        exercise(BaSTIv2.BaSTIv2Library(ts, (feh = [-1.05, -0.9], α_fe = 0.0, canonical = false, diffusion = true, yp = 0.247, η = 0.3)))
    end
    return nothing
end

@compile_workload _precompile_workload()

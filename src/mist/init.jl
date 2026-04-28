
function parse_mist_header(data::AbstractString)
    # Read full file into a string, split by newline
    data = split(data, '\n')
    Y, Z, FeH, aFe, vvcrit = parse.(Float64, filter(!isempty, split(data[5], isspace)[begin+1:end]))
    M, n_pts, n_eep, n_col, phase, type = filter(!isempty, split(data[8], isspace))[begin+1:end]
    M = parse(Float64, M)
    n_pts = parse(Int, n_pts)
    n_eep = parse(Int, n_eep)
    n_col = parse(Int, n_col)
    table_header = filter(!isempty, String.(split(data[12], isspace)))[begin+1:end]
    # MIST readme document says surface metal content should be :log_surf_cell_z,
    # but some files still have the older :log_surf_z keyword. To enforce homogeneity, we will
    # change all occurences of :log_surf_z to :log_surf_cell_z
    idx = findfirst(Base.Fix1(==, "log_surf_z"), table_header)
    if !isnothing(idx)
        table_header[idx] = "log_surf_cell_z"
    end
    return (M_ini = M, Y_ini = Y, Z_ini = Z, feh = FeH, aFe = aFe,
            n_pts = n_pts, n_eep = n_eep, n_col = n_col, phase = phase, type = type,
            colnames=table_header)
end

"""
    mist_feh(fname::AbstractString)
Returns [Fe/H] from the name of MIST track directory.
```jldoctest
julia> StellarTracks.MIST.mist_feh("MIST_v1.2_feh_m0.50_afe_p0.0_vvcrit0.0_EEPS")
-0.5
```
""" 
function mist_feh(fname::AbstractString)
    fname = basename(fname) # Get file name part of path
    feh = parse(Float64, fname[16:19])
    if fname[15] == 'm'
        feh *= -1
    end
    return feh
end
"""
    mist_feh_v2(fname::AbstractString)
Returns [Fe/H] from the name of a MIST v2.5 track archive or directory.
The v2.5 naming convention encodes [Fe/H] as a sign character followed by a
three-digit integer × 100 (e.g., `m400` → -4.00, `p050` → +0.50).
```jldoctest
julia> StellarTracks.MIST.mist_feh_v2("MIST_v2.5_feh_m400_afe_m2_vvcrit0.4_EEPS.txz")
-4.0

julia> StellarTracks.MIST.mist_feh_v2("MIST_v2.5_feh_p050_afe_p0_vvcrit0.0_EEPS.txz")
0.5
```
"""
function mist_feh_v2(fname::AbstractString)
    fname = basename(fname)
    feh = parse(Int, fname[16:18]) / 100
    if fname[15] == 'm'
        feh *= -1
    end
    return feh
end
"""
    mist_afe_v2(fname::AbstractString)
Returns [α/Fe] from the name of a MIST v2.5 track archive or directory.
The v2.5 naming convention encodes [α/Fe] as a sign character followed by a
single integer digit × 10 (e.g., `m2` → -0.2, `p0` → 0.0, `p6` → +0.6).
```jldoctest
julia> StellarTracks.MIST.mist_afe_v2("MIST_v2.5_feh_m400_afe_m2_vvcrit0.4_EEPS.txz")
-0.2

julia> StellarTracks.MIST.mist_afe_v2("MIST_v2.5_feh_p000_afe_p6_vvcrit0.0_EEPS.txz")
0.6
```
"""
function mist_afe_v2(fname::AbstractString)
    fname = basename(fname)
    afe = parse(Int, fname[25:25]) / 10
    if fname[24] == 'm'
        afe *= -1
    end
    return afe
end

"""
    _afe_tag(afe::Number)
Format an [α/Fe] value as the sign+magnitude string used in MIST v2.5 DataDep names,
e.g. `-0.2` → `"m0.2"`, `0.0` → `"p0.0"`, `0.4` → `"p0.4"`.
"""
_afe_tag(afe::Number) = string(afe < 0 ? 'm' : 'p') * string(abs(afe))

"""
    _afe_file_tag(afe::Number)
Format an [α/Fe] value as the single-digit integer tag used in MIST v2.5 archive
filenames, e.g. `-0.2` → `"m2"`, `0.0` → `"p0"`, `0.6` → `"p6"`.
"""
_afe_file_tag(afe::Number) = string(afe < 0 ? 'm' : 'p') * string(round(Int, abs(afe) * 10))

"""
    mist_mass(fname::AbstractString)
Returns initial stellar mass from the name of a single MIST track.
```jldoctest
julia> StellarTracks.MIST.mist_mass("22500M.track.eep")
225.0
```
"""
mist_mass(fname::AbstractString) = parse(Float64, split(basename(fname), '.')[1][begin:end-1]) / 100

"""
    read_mist_track(data::AbstractString; select = nothing)
Given the content of a MIST ".eep" track file as a `String`, parse into
a `TypedTables.Table`. `select` is passed through to `CSV.read` and can be
used to select only certain columns.
"""
function read_mist_track(data::AbstractString; select = nothing)
    header = parse_mist_header(data)
    # This works but is kind of slow; ~4--5 ms compared to the 100 μs to do read(filename, String)
    tdata = Table(CSV.File(IOBuffer(data); comment="#", delim=' ', ignorerepeated=true,
                           header=header.colnames, select=select, ntasks=1))
end
"""
    track_table(filename::AbstractString; select = select_columns))
Given the path to a MIST ".eep" track, read the file as a `TypedTables.Table`.
`select` is passed through to `CSV.read` and can be used to select only certain columns.
"""
function track_table(filename::AbstractString; select = SVector(select_columns))
    result = read_mist_track(read(filename, String); select = select)
    if :Mbol in select
        @argcheck :log_L in select
        # log_L is in solar luminosities, and MIST assumes solar Mbol = 4.74
        mbol = Mbol.(result.log_L, 4.74)
        result = Table(result, Mbol = mbol)
    end
    return result
end
# MIST .eep files have initial whitespace that causes readdlm to fail, see
# https://github.com/JuliaData/DelimitedFiles.jl/issues/1
# track_matrix(filename::AbstractString) = DelimitedFiles.readdlm(filename, ' ', track_type, '\n'; skipstart=13)

function custom_unpack(fname::AbstractString)
    # println(fname)
    # ~/.julia/scratchspaces/124859b0-ceae-595e-8997-d05f6a7a8dfe/datadeps
    # the directory in scratchspaces is the UUID of datadeps
    fpath = dirname(fname)
    fbasename = splitext(basename(fname))[1]
    @info "Unpacking $fbasename"
    out_dir = joinpath(fpath, fbasename)
    if isdir(out_dir)
        rm(out_dir; force=true, recursive=true)
    end
    # unpack_txz is imported from BolometricCorrections.MIST
    # It extracts the .txz, which is an Xz compressed tarball, into the out_dir
    unpack_txz(fname, out_dir)
    
    # Now repack the .track.eep files, which are basic CSV, into JLD2
    feh = string(mist_feh(fname))
    save_dir = joinpath(fpath, feh)
    if isdir(save_dir)
        rm(save_dir; force=true, recursive=true)
    end
    mkdir(save_dir)
    files = filter(Base.Fix1(occursin, ".eep"), readdir(joinpath(out_dir, fbasename); join=true))
    for track in files
        data = read(track, String)
        header = parse_mist_header(data)
        tdata = read_mist_track(data; select = SVector(select_columns))
        JLD2.save_object(joinpath(save_dir, splitext(basename(track))[1]) * ".jld2", tdata)
        # This takes up much more space because of multiple objects I guess?
        # JLD2.jldsave(joinpath(save_dir, splitext(basename(track))[1]) * ".jld2";
        #              table = tdata, M_ini = header.M_ini, Y_ini = header.Y_ini, Z_ini = header.Z_ini,
        #              aFe = header.aFe, n_pts = header.n_pts, n_eep = header.n_eep, n_col = header.n_col,
        #              phase = header.phase, type = header.type)
        # tdata = Table(tdata, m_ini = fill(header.M_ini, length(tdata)))
        # push!(alldata, tdata)
    end
    # Remove temporary directory where .track.eep files were extracted
    rm(out_dir; force=true, recursive=true)
    rm(fname) # Remove original .txz file
end


function custom_unpack_v2(fname::AbstractString)
    fpath = dirname(fname)
    fbasename = splitext(basename(fname))[1]  # e.g. "MIST_v2.5_feh_m400_afe_m2_vvcrit0.4_EEPS"
    @info "Unpacking $fbasename"
    out_dir = joinpath(fpath, fbasename)
    if isdir(out_dir)
        rm(out_dir; force=true, recursive=true)
    end
    # unpack_txz is imported from BolometricCorrections.MIST
    unpack_txz(fname, out_dir)

    # In v2.5 the inner archive directory is the fbasename with the "MIST_vX.Y_" prefix
    # and "_EEPS" suffix stripped (e.g. "feh_m400_afe_m2_vvcrit0.4"), and the .eep files
    # live in a further "eeps/" subdirectory within it.
    inner_dir = replace(replace(fbasename, r"^MIST_v\d+\.\d+_" => ""), r"_EEPS$" => "")
    eep_dir = joinpath(out_dir, inner_dir, "eeps")

    feh = string(mist_feh_v2(fname))
    save_dir = joinpath(fpath, feh)
    if isdir(save_dir)
        rm(save_dir; force=true, recursive=true)
    end
    mkdir(save_dir)
    files = filter(Base.Fix1(occursin, ".eep"), readdir(eep_dir; join=true))
    for track in files
        data = read(track, String)
        tdata = read_mist_track(data; select = SVector(select_columns))
        JLD2.save_object(joinpath(save_dir, splitext(basename(track))[1]) * ".jld2", tdata)
    end
    rm(out_dir; force=true, recursive=true)
    rm(fname)
end

function __init__()
    v1_prefix = "https://mist.science/data/tarballs_v1.2/MIST_v1.2_feh_"
    norot_feh_tags = ["m4.00", "m3.50", "m3.00", "m2.50", "m2.00", "m1.75", "m1.50",
                      "m1.25", "m1.00", "m0.75", "m0.50", "m0.25", "p0.00", "p0.25", "p0.50"]
    norot_dl_links = [v1_prefix * s * "_afe_p0.0_vvcrit0.0_EEPS.txz" for s in norot_feh_tags]

    # Register non-rotating models
    register(DataDep("MISTv1.2_vvcrit0.0",
                     """MIST v1.2 stellar evolutionary tracks without rotation \
                     (v/vcrit=0.0). All available evolutionary tracks \
                     from -4 ≤ [Fe/H] ≤ +0.5 dex will be downloaded.""",
                     norot_dl_links,
                     "5871d51838f238eacc3ceaf02a75e82188eee51c5cf92b5decd0d94377ff692c";
                     post_fetch_method = custom_unpack))

    # Register rotating models
    rot_dl_links = [v1_prefix * s * "_afe_p0.0_vvcrit0.4_EEPS.txz" for s in norot_feh_tags]
    register(DataDep("MISTv1.2_vvcrit0.4",
                     """MIST v1.2 stellar evolutionary tracks with rotation \
                     (v/vcrit=0.4). All available evolutionary tracks \
                     from -4 ≤ [Fe/H] ≤ +0.5 dex will be downloaded.""",
                     rot_dl_links,
                     "94a657226ecb08026df7a205551878559962e191c97740bba045eea3eddd960b";
                     post_fetch_method = custom_unpack))

    # Register MIST v2.5 tracks — one DataDep per (vvcrit, [α/Fe]) combination.
    # Each DataDep downloads all 17 [Fe/H] archives for that combination.
    v2_prefix = "https://mist.science/data/tarballs_v2.5/eeps/MIST_v2.5_feh_"
    v2_feh_tags = ["m400", "m350", "m300", "m275", "m250", "m225", "m200", "m175", "m150",
                   "m125", "m100", "m075", "m050", "m025", "p000", "p025", "p050"]
    # Dict of hashes for the datadeps, keyed by vvcrit and [α/Fe] values as strings (e.g. "0.0", "-0.2", "0.4")
    v2_hashes = Dict("0.0" =>
        Dict("0.0" => "d703498e3c2d544a0dad9ff5663fca0a723ce99f21989061e0ad7b80e8b93bf9",),
        )
    for vvcrit_str in ("0.0", "0.4")
        for afe in afe_grid_v2
            afe_ft  = _afe_file_tag(afe)   # e.g. "m2", "p0", "p2"
            afe_nt  = _afe_tag(afe)        # e.g. "m0.2", "p0.0", "p0.2"
            dl_links = [v2_prefix * feh_tag * "_afe_$(afe_ft)_vvcrit$(vvcrit_str)_EEPS.txz"
                        for feh_tag in v2_feh_tags]
            hash = if vvcrit_str in keys(v2_hashes)
                if string(afe) in keys(v2_hashes[vvcrit_str])
                    v2_hashes[vvcrit_str][string(afe)]
                else
                    "" # No hash available for this combination of vvcrit and [α/Fe] yet
                end
            else
                ""
            end
            println(hash)
            register(DataDep("MISTv2.5_vvcrit$(vvcrit_str)_afe_$(afe_nt)",
                             """MIST v2.5 stellar evolutionary tracks with v/vcrit=$(vvcrit_str) \
                             and [α/Fe]=$(afe). All 17 available [Fe/H] values \
                             (-4 ≤ [Fe/H] ≤ +0.5 dex) will be downloaded.""",
                             dl_links, hash;
                             post_fetch_method = custom_unpack_v2))
        end
    end
end

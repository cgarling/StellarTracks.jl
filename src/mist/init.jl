
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


function __init__()
    norot_feh_tags = ["m4.00", "m3.50", "m3.00", "m2.50", "m2.00", "m1.75", "m1.50",
                      "m1.25", "m1.00", "m0.75", "m0.50", "m0.25", "p0.00", "p0.25", "p0.50"]
    norot_dl_links = ["https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_feh_" * s * "_afe_p0.0_vvcrit0.0_EEPS.txz" for s in norot_feh_tags]

    # Register non-rotating models
    register(DataDep("MISTv1.2_vvcrit0.0",
                     """MIST v1.2 stellar evolutionary tracks without rotation \
                     (v/vcrit=0.0). All available evolutionary tracks \
                     from -4 ≤ [Fe/H] ≤ +0.5 dex will be downloaded.""",
                     norot_dl_links,
                     "5871d51838f238eacc3ceaf02a75e82188eee51c5cf92b5decd0d94377ff692c";
                     post_fetch_method = custom_unpack))

    # Register rotating models
    rot_dl_links = ["https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_feh_" * s * "_afe_p0.0_vvcrit0.4_EEPS.txz" for s in norot_feh_tags]
    register(DataDep("MISTv1.2_vvcrit0.4",
                     """MIST v1.2 stellar evolutionary tracks with rotation \
                     (v/vcrit=0.4). All available evolutionary tracks \
                     from -4 ≤ [Fe/H] ≤ +0.5 dex will be downloaded.""",
                     rot_dl_links,
                     "94a657226ecb08026df7a205551878559962e191c97740bba045eea3eddd960b";
                     post_fetch_method = custom_unpack))
end

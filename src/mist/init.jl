# _parse_mist_linefive(line) = parse.(Float64, filter(!isempty, split(line, isspace)[begin+1:end]))
# function parse_mist_header(fname::AbstractString) # 264.563 μs per file
#     Y, Z, FeH, aFe, vvcrit = 0.0, 0.0, 0.0, 0.0, 0.0
#     for (linenum, line) in enumerate(eachline(fname))
#         if linenum == 5
#             Y, Z, FeH, aFe, vvcrit = _parse_mist_linefive(line)
#             # header = filter(!isempty, split(split(line, string(last(_mist_dependents)))[2], ' '))
#             # return vcat(SVector(string.(_mist_dependents)), convert(Vector{String}, header))
#         end
#     end
#     return (Y_ini = Y, Z_ini = Z, feh = FeH, aFe = aFe)
# end
# This reads the file; new function will take in file content as a string instead
# function parse_mist_header(fname::AbstractString) # 193.867 μs per file
#     # Read full file into a string, split by newline
#     data = split(read(fname, String), '\n')
#     Y, Z, FeH, aFe, vvcrit = parse.(Float64, filter(!isempty, split(data[5], isspace)[begin+1:end]))
#     M, n_pts, n_eep, n_col, phase, type = filter(!isempty, split(data[8], isspace))[begin+1:end]
#     M = parse(Float64, M)
#     n_pts = parse(Int, n_pts)
#     n_eep = parse(Int, n_eep)
#     n_col = parse(Int, n_col)
#     return (M_ini = M, Y_ini = Y, Z_ini = Z, feh = FeH, aFe = aFe,
#             n_pts = n_pts, n_eep = n_eep, n_col = n_col, phase = phase, type = type)
# end
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
    return (M_ini = M, Y_ini = Y, Z_ini = Z, feh = FeH, aFe = aFe,
            n_pts = n_pts, n_eep = n_eep, n_col = n_col, phase = phase, type = type,
            colnames=table_header)
end

# Get [Fe/H] from name of MIST track directory
function mist_feh(fname::AbstractString)
    fname = basename(fname) # Get file name part of path
    feh = parse(Float64, fname[16:19])
    if fname[15] == 'm'
        feh *= -1
    end
    return feh
end

function custom_unpack(fname::AbstractString)
    # println(fname)
    # ~/.julia/scratchspaces/124859b0-ceae-595e-8997-d05f6a7a8dfe/datadeps
    # the directory in scratchspaces is the UUID of datadeps
    fpath = dirname(fname)
    fbasename = splitext(basename(fname))[1]
    out_dir = joinpath(fpath, fbasename)
    if isdir(out_dir)
        rm(out_dir; force=true, recursive=true)
    end
    # unpack_txz is imported from BolometricCorrections.MIST
    # It extracts the .txz, which is an Xz compressed tarball, into the out_dir
    unpack_txz(fname, out_dir)
end

"""
    read_mist_track(data::AbstractString, select = nothing)
Given the content of a MIST ".eep" track file as a `String`, parse into
a `TypedTables.Table`. `select` is passed through to `CSV.read` and can be
used to select only certain columns.
"""
function read_mist_track(data::AbstractString, select = nothing)
    header = parse_mist_header(data)
    # This works but is kind of slow; ~4--5 ms compared to the 100 μs to do read(filename, String)
    tdata = Table(CSV.File(IOBuffer(data); comment="#", delim=' ', ignorerepeated=true,
                           header=header.colnames, select = select, ntasks=1))
end
"""
    track_table(filename::AbstractString)
Given the path to a MIST ".eep" track, read the file as a `TypedTables.Table`.
`select` is passed through to `CSV.read` and can be used to select only certain columns.
"""
function track_table(filename::AbstractString, select = select_columns)
    result = read_mist_track(read(filename, String), select)
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
# track_matrix(filename::AbstractString) = readdlm(filename, ' ', track_type, '\n'; skipstart=13)

function join_tracks(dirname::AbstractString)
    dirs = filter(Base.Fix1(!occursin, ".txz"), readdir(dirname))
    # if length(dirs) != 15
    #     error("Number of track directories not equal to expected 15; \
    # probable unpacking problem.")
    # end
    alldata = []
    for dir in dirs
        # Get the list of .eep track files
        files = filter(Base.Fix1(occursin, ".eep"), readdir(joinpath(dir, dir); join=true))
        for track in files
            data = read(track, String)
            # data = split(rawdata, '\n')
            header = parse_mist_header(data)
            return header
            tdata = Table(CSV.File(IOBuffer(data); comment="#", delim=' ', ignorerepeated=true,
                                                   header=header.colnames, select=select_columns))
            push!(alldata, tdata)
        end
    end
    return alldata
    # Process into one file per unique metallicity, like we did in the PARSEC case
    # a=read("00098M.track.eep", String)
    # Table(CSV.File(IOBuffer(a)))
end


function __init__()
    # dl_links = [
    # "https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_feh_m4.00_afe_p0.0_vvcrit0.4_EEPS.txz",
    # "https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_feh_m3.50_afe_p0.0_vvcrit0.4_EEPS.txz",
    # "https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_feh_m3.00_afe_p0.0_vvcrit0.4_EEPS.txz",
    # "https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_feh_m2.50_afe_p0.0_vvcrit0.4_EEPS.txz",
    # "https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_feh_m2.00_afe_p0.0_vvcrit0.4_EEPS.txz",
    # "https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_feh_m1.75_afe_p0.0_vvcrit0.4_EEPS.txz",
    # "https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_feh_m1.50_afe_p0.0_vvcrit0.4_EEPS.txz",
    # "https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_feh_m1.25_afe_p0.0_vvcrit0.4_EEPS.txz",
    # "https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_feh_m1.00_afe_p0.0_vvcrit0.4_EEPS.txz",
    norot_feh_tags = ["m4.00", "m3.50", "m3.00", "m2.50", "m2.00", "m1.75", "m1.50",
                      "m1.25", "m1.00", "m0.75", "m0.50", "m0.25", "p0.00", "p0.25", "p0.50"]
    norot_dl_links = ["https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_feh_" * s * "_afe_p0.0_vvcrit0.0_EEPS.txz" for s in norot_feh_tags]
    
    register(DataDep("MISTv1.2_vvcrit0.0",
                     """MIST v1.2 stellar evolutionary tracks without rotation \
                     (v/vcrit=0.0). All available evolutionary tracks \
                     from -4 ≤ [Fe/H] ≤ +0.5 dex will be downloaded.""",
                     norot_dl_links,
                     "5871d51838f238eacc3ceaf02a75e82188eee51c5cf92b5decd0d94377ff692c"))

end

# First, utilities for parsing the Rosenfield 2016 track files

function file_properties(filename::AbstractString) # extract properties of original track files
    file = basename(filename) # split(filename, "/")[end]
    Z = parse(track_type, split( split(file, "Z")[2], "Y" )[1])
    Y = parse(track_type, split( split(file, "Y")[2], "OUT" )[1])
    M_HB = split( split(file, "M")[2], ".dat" )[1]
    if occursin("HB", M_HB)
        HB = true
        M = parse(track_type, split(M_HB, ".HB")[1])
    else
        HB = false
        M = parse(track_type, M_HB)
    end
    return (Z=Z, Y=Y, M=M, HB=HB)
end

"""
    dir_properties(dir::AbstractString, checkdir::Bool=true)
Parse `Z` and `Y` from name of track directory (e.g., ".../tracks/Z0.0001Y0.249/") with or without trailing '/'.
If `checkdir` is `true`, this function will check that the provided directory `dir` exists on the file system. """
function dir_properties(dir::AbstractString, checkdir::Bool=true) # utility
    if checkdir
        @argcheck isdir(dir) ArgumentError("`checkdir` argument to `dir_properties` is `true` and file system directory `dir` does not exist.")
    end
    if dir[end] == '/'
        dir = split(dir, '/')[end-1]
    else
        dir = split(dir, '/')[end]
    end
    Z = parse(track_type, split(split(dir, 'Z')[end], 'Y')[1])
    Y = parse(track_type, split(dir, 'Y')[end])
    return (Z=Z, Y=Y)
end

# CSV actually seems faster, but this method returns a matrix rather than Table like CSV.
"""
    track_matrix(filename::AbstractString)
Given the path to a PARSEC v1.2S ".dat" track, read the file as a `Matrix`.
"""
track_matrix(filename::AbstractString) = readdlm(filename, ' ', track_type, '\n'; skipstart=1)

"""
    track_table(filename::AbstractString, select = StellarTracks.PARSEC.select_columns)
Given the path to a PARSEC v1.2S ".dat" track, read the file as a `TypedTables.Table`.
`select` is passed through to `CSV.read` and can be used to select only certain columns.
"""
function track_table(filename::AbstractString, select = select_columns)
    return CSV.read(filename, Table;
                    skipto = 2, header = Vector(track_header), types = track_type,
                    select = select) # utility
end

#########################################################################

# Extract gzipped tar (.tar.gz)
function ungzip(fname::AbstractString, dir::AbstractString)
    @argcheck isfile(fname)
    gz = open(fname)
    tar = GzipDecompressorStream(gz)
    # Tar.extract(tar, joinpath(dirname(fname), dir))
    Tar.extract(tar, String(dir))
    close(tar)
end

# function custom_unpack(fpath::AbstractString=datadep"PARSECv1.2S")
function custom_unpack(fname::AbstractString=joinpath(datadep"PARSECv1.2S", "releasev2.tar.gz"))
    @info "Unpacking PARSEC v1.2S stellar tracks"
    fpath = dirname(fname)
    # fname = joinpath(fpath, "releasev2.tar.gz")
    # Extract gzipped tar; delete existing directory if it exists
    out_dir = joinpath(fpath, "releasev2")
    if isdir(out_dir)
        rm(out_dir; force=true, recursive=true)
    end
    ungzip(fname, out_dir)
    # Uncompress containing directories
    files = readdir(out_dir; join=true)
    for fi in files
        tmpdir = split(fi, ".tar.gz")[1]
        ungzip(fi, tmpdir)
        rm(fi)
    end

    @info "Reorganizing PARSEC v1.2S stellar tracks"
    # Process individual files into joined files
    @showprogress for dir in readdir(out_dir; join=true)
        println("First data directory:\n")
        dir |> println
        # Z, Y = dir_properties(dir)
        Z = parse(Float64, split(split(dir, 'Z')[end], 'Y')[1])
        Y = parse(Float64, split(dir, 'Y')[end])
        # dirstem = split(dir, "/")[end] # Something like Z0.0001Y0.249
        dirstem = splitpath(dir)[end] # Something like Z0.0001Y0.249
        dir = joinpath(dir, dirstem)
        println("Joined data directory:\n")
        dir |> println
        files = readdir(dir; join=true) # List of files in directory
        println("All files\n")
        files |> display
        # Remove ._ files
        files = filter(Base.Fix1(!occursin, "._"), files)
        # Figure out which files are HB only and which contain main sequence -> trgb
        hb_idx = occursin.(".HB.", files) # BitVector
        hb_files = files[hb_idx]
        ms_files = files[broadcast(!, hb_idx)]

        println("MS files\n")
        ms_files |> display
        println("HB files\n")
        hb_files |> display

        # Load all MS files and process into a single data table for more efficient access        
        ms_props = Table(file_properties(file) for file in ms_files)
        # Sort files according to initial mass
        idxs = sortperm(ms_props.M)
        ms_props = Table(ms_props[idxs])
        ms_files = ms_files[idxs]
        # Construct output containers
        data = Matrix{track_type}(undef, 0, 6) # length(keepcols)) # Doing away with keepcols
        masses = Vector{track_type}(undef, 0) # Keep track of initial masses
        eep_vec = Vector{Int}(undef, 0)       # Keep track of EEP points
        for (i, file) in enumerate(ms_files)
            mass = file_properties(file).M
            tmpdata = track_matrix(file)
            # tmpdata = view(track_matrix(file), :, keepcols) # Doing away with keepcols
            data = vcat(data, tmpdata)
            masses = vcat(masses, fill(mass, size(tmpdata, 1)))
            eep_vec = vcat(eep_vec, 1:size(tmpdata, 1))
        end
        # data = Table(NamedTuple{track_header_symbols[keepcols]}(row) for row in eachrow(data))
        data = Table(NamedTuple{track_header_symbols}(row) for row in eachrow(data))
        data = Table(data, m_ini = masses, eep = eep_vec)

        # To use the HB files, we need to identify the age for each track at the TRGB
        # as the logAge entries in the .HB files start from the ZAHB, and the .HB files
        # contain tracks with stellar masses that are not computed in the MS files.
        # Therefore need to interpolate x=m_ini against y=logAge at the TRGB.
        hbdata = Matrix{track_type}(undef, 0, 6) # length(keepcols)) # Doing away with keepcols
        tmpdata = data[data.eep .== (eep_idxs.RG_TIP - 1)]
        tmpmrange = extrema(tmpdata.m_ini)
        trgb_interp = AkimaInterpolation(tmpdata.logAge, tmpdata.m_ini)
        la_col = findfirst(track_header .== "logAge") # Find column that contains logAge

        # Construct output containers
        trgb_logAge = Vector{track_type}(undef, 0)
        hb_masses = similar(trgb_logAge)   # Keep track of initial masses
        hb_eep_vec = Vector{Int}(undef, 0) # Keep track of EEP points
        for (i, file) in enumerate(hb_files)
            fprops = file_properties(file)
            mass = fprops.M
            if tmpmrange[1] <= mass <= tmpmrange[2] # If mass within interpolation range,
                trgb_la = trgb_interp(mass)
                push!(trgb_logAge, trgb_la)
                # Load hb track data
                # tmphbdata = view(track_matrix(file), :, keepcols) # Doing away with keepcols
                tmphbdata = track_matrix(file)
                # Add ZAHB logAge in file to TRGB logAge calculated above
                # tmphbdata = Table(tmphbdata, logAge = log10.(exp10.(tmphbdata.logAge) .+ exp10(trgb_la)))
                @views tmphbdata[:, la_col] .= log10.(exp10.(tmphbdata[:, la_col]) .+ exp10(trgb_la))
                hbdata = vcat(hbdata, tmphbdata)
                hb_masses = vcat(hb_masses, fill(mass, size(tmphbdata, 1)))
                hb_eep_vec = vcat(hb_eep_vec, (0:size(tmphbdata, 1)-1) .+ eep_idxs.RG_TIP) # Offset EEP appropriately
            end
        end
        # hbdata = Table(NamedTuple{track_header_symbols[keepcols]}(row) for row in eachrow(hbdata))
        hbdata = Table(NamedTuple{track_header_symbols}(row) for row in eachrow(hbdata))
        hbdata = Table(hbdata, m_ini = hb_masses, eep = hb_eep_vec)
        data = vcat(data, hbdata)
        # 1/3 size of uncompressed CSV, 2x read time, 80 ms per file -- totally fine
        # CSV.write(joinpath(fpath, dirstem)*".gz", data; compress=true)
        # JLD2 loading takes 1 ms per file -- much faster than CSV ...
        JLD2.save_object(joinpath(fpath, dirstem * ".jld2"), data)
        # ; compress=true) halves size, 25x load time
        rm(dir; recursive=true)
    end
    @info "Cleaning up data directory"
    rm(fname) # remove releasev2.tar.gz
    rm(out_dir; recursive=true)
end

function __init__()
    register(DataDep("PARSECv1.2S",
"""
PARSEC V1.2S stellar evolutionary tracks with
equivalent evolutionary points calculated by Phil Rosenfield.

This is specifically the v2.0 release from the repository
https://github.com/philrosenfield/padova_tracks
with permanent DOI
10.5281/zenodo.48123.

See the above links for information on how to cite these data.

This data product is licensed under the following terms:

The MIT License (MIT)

Copyright (c) 2015 Phil Rosenfield

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. """,
                     "https://github.com/philrosenfield/padova_tracks/releases/download/v2.0/releasev2.tar.gz",
                     "fc104dd040c3a8b4dc9d71d3313ef35efa5eecbb93a63884b2a910532db07e5e";
                     post_fetch_method = custom_unpack))
                     # ))
end

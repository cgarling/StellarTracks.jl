module MIST

# imports from parent module
using ..StellarTracks: AbstractTrack, AbstractTrackSet, AbstractTrackLibrary, uniqueidx, Mbol
import ..StellarTracks: mass, post_rgb, isochrone
import ..StellarTracks: X, X_phot, Y, Y_phot, Z, Z_phot, MH, chemistry

# Imports for data reading / processing
import CSV
using DataDeps: register, DataDep, @datadep_str, unpack
using DataInterpolations: AbstractInterpolation, CubicSpline, CubicHermiteSpline, PCHIPInterpolation
# using DelimitedFiles: readdlm
using Glob: glob
import JLD2 # for saving files in binary format
using ProgressMeter: @showprogress
using TypedTables: Table

# We will simply reuse the MISTChemistry defined in BolometricCorrections
using BolometricCorrections.MIST: MISTChemistry, unpack_txz

# Imports for core module code
using ArgCheck: @argcheck
using StaticArrays: SVector

# const eep_idxs = (PMS_BEG = 1,      # beginning of PMS; log(T_c) = 5
#                   MS_BEG = 202,     # beginning of MS; H-burning L > 99.9% of total L
#                   IAMS = 353,       # intermediate age MS; X_c = 0.3
#                   MS_TO = 454,      # terminal age MS; X_c = 1e-12
#                   RG_TIP = 605,     # max(L) or min(Teff) after core H-burning, before core-He burning
#                   HE_BEG = 631,     # beginning core He fusion to end core He fusion
#                   END_CHEB = 707,   # end of core He burning; Y_c = 1e-4
#                   TPAGB_BEG = 808,  # For low mass stars, this is the beginning of TP-AGB;
#                                     # Y_c < 1e-6, H and He shells comparable mass;
#                                     # for high mass stars, this is the end of core carbon burning
#                   END_TPAGB = 1409, # For stars that will become white dwarfs, end of TP-AGB phase
#                   WD = 1710)        # White dwarfs
# # Indices where the different phases begin

# const eep_lengths = NamedTuple{keys(eep_idxs)[begin:end-1]}(eep_idxs[i+1] - eep_idxs[i] for i in eachindex(values(eep_idxs))[begin:end-1])

const eep_lengths = (PMS_BEG = 201,   # beginning of PMS; log(T_c) = 5
                     MS_BEG = 151,    # beginning of MS; H-burning L > 99.9% of total L
                     IAMS = 101,      # intermediate age MS; X_c = 0.3
                     MS_TO = 151,     # terminal age MS; X_c = 1e-12
                     RG_TIP = 26,     # max(L) or min(Teff) after core H-burning, before core-He burning
                     HE_BEG = 76,     # beginning core He fusion to end core He fusion
                     END_CHEB = 101,  # end of core He burning; Y_c = 1e-4
                     TPAGB_BEG = 601, # For low mass stars, this is the beginning of TP-AGB;
                                      # Y_c < 1e-6, H and He shells comparable mass;
                                      # for high mass stars, this is the end of core carbon burning
                     END_TPAGB = 301) # For stars that will become white dwarfs, end of TP-AGB phase
                                      # White dwarfs at 1710, final EEP point
# Indices where the different phases begin
const eep_idxs = NamedTuple{keys(eep_lengths)}((1, (cumsum(values(eep_lengths)[begin:end-1]) .+ 1)...))
# Which columns to actually keep after reading track file; used in Track and track_table
const select_columns = SVector(:star_age, :log_L, :log_Teff, :log_g, :log_surf_cell_z) # :star_mass,
const track_type = Float64 # MIST track files have Float64 precision
const feh_grid = SVector(-4.0, -3.5, -3.0, -2.75, -2.5, -2.25, -2.0, -1.75, -1.5, -1.25, -1.0,
                         -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75)


##########################################################################

# Data download, organization, etc.
include("init.jl")

##########################################################################

""" `MISTTrack` implements the [`AbstractTrack`](@ref StellarTracks.AbstractTrack)
interface for the MIST stellar evolution library. """
struct MISTTrack{A,B,C} <: AbstractTrack
    filename::String
    data::Table{A}
    itp::B
    properties::C
end
function MISTTrack(feh::Number, mass::Number, vvcrit::Number=0)
    props = (M = mass, feh = feh)
    # Validate feh
    @argcheck feh in feh_grid
    feh = string(feh_grid[searchsortedfirst(feh_grid, feh)])
    # Validate vvcrit
    if vvcrit == 0
        vvcrit = "0.0"
    elseif vvcrit == 0.4
        vvcrit = "0.4"
    else
        throw(ArgumentError("Invalid vvcrit=$vvcrit argument; valid arguments are 0 or 0.4."))
    end
    dd_path = @datadep_str("MISTv1.2_vvcrit"*vvcrit)
    # Validate mass
    allfiles = readdir(joinpath(dd_path, feh); join=true)
    masses = mist_mass.(allfiles)
    @argcheck mass in masses ArgumentError("Invalid mass=$mass argument; available track masses for [Fe/H]=$feh and vvcrit=$vvcrit are $masses.")
    mass = string(Int(mass * 100))
    mass = repeat("0", 5 - length(mass)) * mass # Pad to length 5
    # Load data file into table
    filename = joinpath(dd_path, feh, mass * "M.track.jld2")
    data = JLD2.load_object(filename)
    # return data
    # Construct interpolator as a function of proper age
    # itp = interpolate(data.star_age,
    #                   [SVector(values(i)[2:end]) for i in data], Gridded(Linear()))
    itp = CubicSpline([SVector(values(i)[2:end]) for i in data], data.star_age)
    return MISTTrack(convert(String, filename), data, itp, props)
end
# Make Track callable with logAge to get logTe, Mbol, and logg as a NamedTuple
function (track::MISTTrack)(logAge::Number)
    # result = track.itp(logAge)
    result = track.itp(exp10(logAge))
    return NamedTuple{Tuple(select_columns)[2:end]}(result)
end
Base.extrema(t::MISTTrack) = log10.(extrema(t.itp.t))
mass(t::MISTTrack) = t.properties.M
chemistry(::MISTTrack) = MISTChemistry()
MH(t::MISTTrack) = t.properties.feh # MH(chemistry(t), Z(t))
Z(t::MISTTrack) = Z(chemistry(t), MH(t)) # t.properties.Z
Y(t::MISTTrack) = Y(chemistry(t), Z(t))  # t.properties.Y
X(t::MISTTrack) = 1 - Y(t) - Z(t)
# Whether there is post-rgb evolution or not is dependent on how many EEPs in track
post_rgb(t::MISTTrack) = length(t.itp.u) > eep_idxs.RG_TIP
Base.eltype(t::MISTTrack) = typeof(t.properties.feh)

##########################################################################

""" `MISTTrackSet` implements the [`AbstractTrackSet`](@ref StellarTracks.AbstractTrackSet)
interface for the MIST stellar evolution library. """
struct MISTTrackSet{A <: AbstractVector{<:Integer},
                    B <: AbstractVector{<:AbstractInterpolation}, # AbstractInterpolation{T}
                    # C <: AbstractVector{<:AbstractVector{<:AbstractInterpolation}},
                    C,
                    D} <: AbstractTrackSet
    eeps::A # EEP points; vector of integers
    AMRs::B # age-mass relations; vector of interpolators, same length as eeps
    interps::C
    properties::D
end
# Given a metallicity and rotation, load the correct models and then call below method
function MISTTrackSet(feh::Number, vvcrit::Number=0) # One table per stellar model
    chem = MISTChemistry()
    zval = Z(chem, feh)
    # Validate feh
    @argcheck feh in feh_grid
    feh = string(feh_grid[searchsortedfirst(feh_grid, feh)])
    # Validate vvcrit
    if vvcrit == 0
        vvcrit = "0.0"
    elseif vvcrit == 0.4
        vvcrit = "0.4"
    else
        throw(ArgumentError("Invalid vvcrit=$vvcrit argument; valid arguments are 0 or 0.4."))
    end
    dd_path = @datadep_str("MISTv1.2_vvcrit"*vvcrit)
    # List stellar track files
    allfiles = readdir(joinpath(dd_path, feh); join=true)
    masses = mist_mass.(allfiles)
    # return [JLD2.load_object(file) for file in allfiles]
    # return MISTTrackSet([JLD2.load_object(file) for file in allfiles],
    #                     masses, parse(track_type, feh))
    data = vcat([begin
                     tmpdata = JLD2.load_object(allfiles[i])
                     Table(tmpdata, m_ini=fill(masses[i], length(tmpdata)),
                           eep = 1:length(tmpdata))
                 end for i in eachindex(allfiles, masses)]...)
    # return data
    return MISTTrackSet(data, parse(track_type, feh))
end
# One table per stellar model
# function MISTTrackSet(data::Vector{<:Table}, masses::Vector{<:Number}, feh::Number)
#     eeps = 1:maximum(length.(data))
#     itp_type = CubicHermiteSpline{Vector{track_type},
#                                   Vector{track_type},
#                                   Vector{track_type},
#                                   Vector{track_type},
#                                   Vector{track_type},
#                                   track_type}
#     amrs = Vector{itp_type}(undef, length(eeps))
#     logte = Vector{itp_type}(undef, length(eeps))
#     logl = similar(logte)
#     logg = similar(logte)
#     logsurfz = similar(logte)

#     Threads.@threads for i in eachindex(eeps)
#         tmpdata = Table(data[i]
#     end
# end
function MISTTrackSet(data::Table, feh::Number)
    # eeps = 1:maximum(length.(data))
    eeps = sort(unique(data.eep))
    itp_type = CubicHermiteSpline{Vector{track_type},
                                  Vector{track_type},
                                  Vector{track_type},
                                  Vector{track_type},
                                  Vector{track_type},
                                  track_type}
    amrs = Vector{itp_type}(undef, length(eeps))
    logte = Vector{itp_type}(undef, length(eeps))
    logl = similar(logte)
    logg = similar(logte)
    logsurfz = similar(logte)

    # Threads.@threads for i in eachindex(eeps)
    for i in eachindex(eeps)
        eep = eeps[i]
        tmpdata = data[findall(Base.Fix1(==, eep), data.eep)] # Performance optimization
        # Sort by age, which will be the independent variable in the AMR interpolation
        idxs = sortperm(tmpdata.star_age)
        tmpdata = tmpdata[idxs]
        # Now keep only entries that define a unique mapping between age and m_ini
        tmpdata = tmpdata[uniqueidx(tmpdata.star_age)]

        # We'll start the interpolation at the most massive model and enforce a monotonic
        # decrease in stellar initial mass with increasing age at fixed EEP
        _, p1 = findmax(tmpdata.m_ini)
        idxs = p1:lastindex(tmpdata)
        goodidxs = diff(tmpdata.m_ini[idxs]) .< 0
        goodidxs = vcat(true, goodidxs) # add true for first element as well
        tmpdata = tmpdata[idxs[goodidxs]]

        # PCHIPInterpolation is a type of CubicHermiteSpline
        amrs[i] = PCHIPInterpolation(tmpdata.m_ini, log10.(tmpdata.star_age))
        # Sort by initial stellar mass for defining the other interpolations
        idxs = sortperm(tmpdata.m_ini)
        tmpdata = tmpdata[idxs]
        # Now interpolate mass against logte, mbol, logg, c_o
        logte[i] = PCHIPInterpolation(tmpdata.log_Teff, tmpdata.m_ini)
        logl[i] = PCHIPInterpolation(tmpdata.log_L, tmpdata.m_ini)
        logg[i] = PCHIPInterpolation(tmpdata.log_g, tmpdata.m_ini)
        logsurfz[i] = PCHIPInterpolation(tmpdata.log_surf_cell_z, tmpdata.m_ini)
    end
    # chem = MISTChemistry()
    # zval = Z(chem, feh)
    # return MISTTrackSet(eeps, amrs, (log_L = logl, log_Teff = logte, log_g = logg, log_surf_cell_z = logsurfz), (Z = zval, Y = Y(chem, zval), masses = unique(data.m_ini)))
    return MISTTrackSet(eeps, amrs, (log_L = logl, log_Teff = logte, log_g = logg, log_surf_cell_z = logsurfz), (feh = feh, masses = unique(data.m_ini)))
end
function (ts::MISTTrackSet)(M::Number) # Interpolation to get a Track with mass M
    error("Not yet implemented.")
end
mass(ts::MISTTrackSet) = ts.properties.masses
chemistry(::MISTTrackSet) = MISTChemistry()
MH(ts::MISTTrackSet) = ts.properties.feh # MH(chemistry(ts), Z(ts))
Z(ts::MISTTrackSet) = Z(chemistry(ts), MH(ts)) # ts.properties.Z
Y(ts::MISTTrackSet) = Y(chemistry(ts), Z(ts)) # ts.properties.Y
X(ts::MISTTrackSet) = 1 - Y(ts) - Z(ts)
Base.eltype(ts::MISTTrackSet) = typeof(ts.properties.feh)
function Base.show(io::IO, mime::MIME"text/plain", ts::MISTTrackSet)
    print(io, "MISTTrackSet with Y=$(Y(ts)), Z=$(Z(ts)), $(length(ts.AMRs)) EEPs and $(length(ts.properties.masses)) initial stellar mass points.")
end

function isochrone(ts::MISTTrackSet, logAge::Number) # 1 ms
    eeps = Vector{Int}(undef, 0)
    track_extrema = extrema(ts.properties.masses)
    interp_masses = Vector{eltype(ts)}(undef, 0)
    logte = similar(interp_masses)
    logl = similar(interp_masses)
    logg = similar(interp_masses)
    logsurfz = similar(interp_masses)
    # data = Vector{eltype(first(ts.AMRs))}(undef, 0)
    for (i, amr) in enumerate(ts.AMRs)
        drange = extrema(amr) # Get valid age range for the EEP
        if logAge >= first(drange) && logAge <= last(drange)
            imass = amr(logAge) # 100 ns
            # imass = amr(age)
            # Due to non-linear interpolation of mass(age), it is possible for the
            # interpolated masses to be outside the range valid for the tracks
            # Check here and do not write if outside range. Also require monotonically
            # increasing initial masses (gets rid of some non-monotonic behavior in the
            # underlying EEP tracks).
            if imass >= first(track_extrema) && imass <= last(track_extrema)
                logli = ts.interps.log_L[i](imass) # 120 ns
                # Enforce monotonically increasing luminosity along the MS (which ends at eep_idxs[4])
                # We can't be sure if last point was overly bright or current point is overly faint,
                # so we delete last point and continue so this point isn't output
                # Maybe not necessary, revisit later
                if length(logl) > 0 && i < eep_idxs[4] && logli < last(logl)
                    li = lastindex(logl)
                    deleteat!(eeps, li)
                    deleteat!(interp_masses, li)
                    deleteat!(logte, li)
                    deleteat!(logl, li)
                    deleteat!(logg, li)
                    deleteat!(logsurfz, li)
                    continue
                end
                push!(eeps, i)
                push!(interp_masses, imass)
                push!(logl, logli)
                push!(logte, ts.interps.log_Teff[i](imass))
                push!(logg, ts.interps.log_g[i](imass))
                push!(logsurfz, ts.interps.log_surf_cell_z[i](imass))
            end
        end
    end
    # return (eep = eeps, m_ini = interp_masses, log_L = logl, log_Teff = logte, log_g = logg, log_surf_cell_z = logsurfz)
    return (eep = eeps, m_ini = interp_masses, logTe = logte, Mbol = Mbol.(logl, 4.74),
            logg = logg, logL = logl, log_surf_cell_z = logsurfz)
end

# export PARSECLibrary, PARSECChemistry, MH_canon, Z_canon # Unique module exports
export MISTTrack, MISTTrackSet # Unique module exports
export mass, X, Y, Z, MH, post_rgb, isochrone # Export generic API methods

end # module

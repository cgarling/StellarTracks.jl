"""StellarTracks.MIST provides access to the MIST stellar tracks."""
module MIST

# imports from parent module
using ..StellarTracks: AbstractTrack, AbstractTrackSet, AbstractTrackLibrary,
                       uniqueidx, Mbol, _generic_trackset_interp
import ..StellarTracks: X, Y, Z, MH, chemistry, mass, post_rgb, isochrone # X_phot, Y_phot, Z_phot,

# Imports for data reading / processing
import CSV
using DataDeps: register, DataDep, @datadep_str, unpack
using DataInterpolations: AbstractInterpolation, CubicSpline, CubicHermiteSpline, PCHIPInterpolation
# # MIST track with [Fe/H] = -4, M=0.9, vvcrit=0.0 has float error in the star_age that results
# in non-increasing as a function of EEP. Since it's ~EPS error, we use Interpolations.deduplicate_knots!
# to correct for this minor error.
using Interpolations: deduplicate_knots!
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

"""Number of secondary EEP points per primary EEP point."""
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
"""1-based indices giving the EEP point at which each EEP phase begins."""
const eep_idxs = NamedTuple{keys(eep_lengths)}((1, (cumsum(values(eep_lengths)[begin:end-1]) .+ 1)...))
# Which columns to actually keep after reading track file; used in Track and track_table
"""Columns to save from the MIST tracks."""
const select_columns = (:star_age, :log_L, :log_Teff, :log_g, :log_surf_cell_z) # :star_mass,
"""Data type to parse the MIST tracks as."""
const track_type = Float64 # MIST track files have Float64 precision
"""Available [Fe/H] values in the MIST stellar track grid."""
const feh_grid = track_type[-4.0, -3.5, -3.0, -2.5, -2.0, -1.75, -1.5, -1.25, -1.0,
                            -0.75, -0.5, -0.25, 0.0, 0.25, 0.5]
"""Initial stellar masses for the stellar tracks in the MIST grid. These are uniform for all [Fe/H] and also for rotating and non-rotating grids."""
const mass_grid = track_type[0.1, 0.15, 0.2, 0.25, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.36, 1.38, 1.4, 1.42, 1.44, 1.46, 1.48, 1.5, 1.52, 1.54, 1.56, 1.58, 1.6, 1.62, 1.64, 1.66, 1.68, 1.7, 1.72, 1.74, 1.76, 1.78, 1.8, 1.82, 1.84, 1.86, 1.88, 1.9, 1.92, 1.94, 1.96, 1.98, 2.0, 2.02, 2.04, 2.06, 2.08, 2.1, 2.12, 2.14, 2.16, 2.18, 2.2, 2.22, 2.24, 2.26, 2.28, 2.3, 2.32, 2.34, 2.36, 2.38, 2.4, 2.42, 2.44, 2.46, 2.48, 2.5, 2.52, 2.54, 2.56, 2.58, 2.6, 2.62, 2.64, 2.66, 2.68, 2.7, 2.72, 2.74, 2.76, 2.78, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0]

"""
    _parse_vvcrit(vvcrit::Number)
Given a MIST vvcrit value, determine if it is valid and return a `String`
corresponding to the proper value.
```jldoctest
julia> using StellarTracks.MIST: _parse_vvcrit

julia> _parse_vvcrit(0.0)
"0.0"

julia> _parse_vvcrit(0.4)
"0.4"

julia> using Test; @test_throws ArgumentError _parse_vvcrit(0.2)
Test Passed
      Thrown: ArgumentError
```
"""
function _parse_vvcrit(vvcrit::Number)
    # Validate vvcrit
    if vvcrit ≈ 0
        return "0.0"
    elseif vvcrit ≈ 0.4
        return "0.4"
    else
        throw(ArgumentError("Invalid vvcrit=$vvcrit argument; valid arguments must be ≈ 0 or 0.4."))
    end
end

"""
    _parse_mass(mass::Number)
Given a MIST stellar initial mass value, determine if it is valid and return a `String` giving the leading part of a track filename specification.

```jldoctest
julia> using StellarTracks.MIST: _parse_mass

julia> _parse_mass(1.0)
"00100"

julia> using Test; @test_throws ArgumentError _parse_mass(1.011)
Test Passed
      Thrown: ArgumentError
```
"""
function _parse_mass(mass::Number)
    mass_idx = findfirst(≈(mass), mass_grid)
    @argcheck !isnothing(mass_idx) ArgumentError("Invalid mass=$mass argument; available track masses are $mass_grid. For track interpolation, use `MISTLibrary`.")
    # mass = string(Int(mass * 100))
    mass = string(round(Int, mass_grid[mass_idx] * 100))
    mass = repeat("0", 5 - length(mass)) * mass # Pad to length 5
    return mass
end

##########################################################################

# Data download, organization, etc.
include("init.jl")

##########################################################################

"""
    MISTTrack(mh::Number, mass::Number, vvcrit::Number=0)
`MISTTrack` implements the [`AbstractTrack`](@ref StellarTracks.AbstractTrack)
interface for the MIST stellar evolution library.
```jldoctest
julia> track = StellarTracks.MIST.MISTTrack(-2, 0.15, 0.0)
MISTTrack with M_ini=0.15, MH=-2, vvcrit=0.0, Z=0.00014899227095992976, Y=0.24922374966753474, X=0.7506272580615054.

julia> track(7.0) # interpolate track at log10(age [yr]) = 7
(log_L = -1.5293719450743, log_Teff = 3.587337261741102, log_g = 4.447603584617846, log_surf_cell_z = -3.8450984758441953)
```
"""
struct MISTTrack{A,B,C} <: AbstractTrack
    data::Table{A}
    itp::B
    properties::C
end
function MISTTrack(data::Table, props)
    # Construct interpolator as a function of proper age
    # itp = interpolate(deduplicate_knots!(data.star_age; move_knots=true),
    #                   [SVector(values(i)[2:end]) for i in data], Gridded(Linear()))
    itp = CubicSpline([SVector(values(i)[2:end]) for i in data],
                      deduplicate_knots!(data.star_age; move_knots=true);
                      cache_parameters=true)
    return MISTTrack(data, itp, props)
end
# Given feh, mass, vvcrit, load the data table and call to function above
function MISTTrack(feh::Number, mass::Number, vvcrit::Number=0)
    props = (M = mass, feh = feh, vvcrit = vvcrit)
    feh_idx = findfirst(≈(feh), feh_grid) # Validate against feh_grid
    @argcheck !isnothing(feh_idx) ArgumentError("Provided `feh` argument $feh to `MISTTrack` is invalid; available metallicities are $feh_grid. For metallicity interpolation, use `MISTLibrary`.")
    feh = string(feh_grid[feh_idx])
    # Validate vvcrit
    vvcrit = _parse_vvcrit(vvcrit)
    # Validate mass
    mass = _parse_mass(mass)
    # Load data file into table
    filename = @datadep_str(joinpath("MISTv1.2_vvcrit"*vvcrit, feh, mass * "M.track.jld2"))
    data = JLD2.load_object(filename)
    return MISTTrack(data, props)
end
# Make Track callable with logAge to get properties as a NamedTuple
Base.keys(::MISTTrack) = select_columns[2:end]
function (track::MISTTrack)(logAge::Number)
    # result = track.itp(logAge)
    result = track.itp(exp10(logAge))
    return NamedTuple{keys(track)}(result)
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
function Base.show(io::IO, mime::MIME"text/plain", t::MISTTrack)
    print(io, "MISTTrack with M_ini=$(mass(t)), MH=$(MH(t)), vvcrit=$(t.properties.vvcrit), Z=$(Z(t)), Y=$(Y(t)), X=$(X(t)).")
end

##########################################################################

"""
    MISTTrackSet(mh::Number, vvcrit::Number=0)
`MISTTrackSet` implements the [`AbstractTrackSet`](@ref StellarTracks.AbstractTrackSet)
interface for the MIST stellar evolution library.
```jldoctest
julia> ts = StellarTracks.MIST.MISTTrackSet(0.0, 0.0)
MISTTrackSet with MH=0.0, vvcrit=0.0, Z=0.0142014201420142, Y=0.2703270327032703, 1710 EEPs and 196 initial stellar mass points.

julia> ts(1.01) # Interpolate track at new initial mass
MISTTrack with M_ini=1.01, MH=0.0, vvcrit=0.0, Z=0.0142014201420142, Y=0.2703270327032703, X=0.7154715471547155.

julia> isochrone(ts, 10.0) isa NamedTuple # Interpolate isochrone at `log10(age [yr]) = 10`
true
```
"""
struct MISTTrackSet{A <: AbstractVector{<:Integer},
                    B <: AbstractVector{<:AbstractInterpolation},
                    C,
                    D} <: AbstractTrackSet
    eeps::A # EEP points; vector of integers
    AMRs::B # age-mass relations; vector of interpolators, same length as eeps
    interps::C
    properties::D
end
# Given a metallicity and rotation, load the correct models and then call below method
function MISTTrackSet(feh::Number, vvcrit::Number=0) # One table per stellar model
    # Validate feh
    feh_idx = findfirst(≈(feh), feh_grid) # Validate against feh_grid
    @argcheck !isnothing(feh_idx) ArgumentError("Provided `feh` argument $feh to `MISTTrackSet` is invalid; available metallicities are $feh_grid. For metallicity interpolation, use `MISTLibrary`.")
    feh = string(feh_grid[feh_idx])
    # Validate vvcrit
    vvcrit = _parse_vvcrit(vvcrit)
    # List stellar track files
    allfiles = [@datadep_str(joinpath("MISTv1.2_vvcrit"*vvcrit, feh, _parse_mass(mass) * "M.track.jld2")) for mass in mass_grid]
    # return [JLD2.load_object(file) for file in allfiles]
    # return MISTTrackSet([JLD2.load_object(file) for file in allfiles],
    #                     masses, parse(track_type, feh))
    data = vcat([begin
                     tmpdata = JLD2.load_object(allfiles[i])
                     Table(tmpdata, m_ini=fill(mass_grid[i], length(tmpdata)),
                           eep = 1:length(tmpdata))
                 end for i in eachindex(allfiles)]...)
    # return data
    return MISTTrackSet(data, parse(track_type, feh), parse(track_type, vvcrit))
end
function MISTTrackSet(data::Table, feh::Number, vvcrit::Number)
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

    Threads.@threads for i in eachindex(eeps)
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
    return MISTTrackSet(eeps, amrs,
                        (log_L = logl, log_Teff = logte, log_g = logg, log_surf_cell_z = logsurfz),
                        (feh = feh, vvcrit = vvcrit, masses = unique(data.m_ini)))
end
function (ts::MISTTrackSet)(M::Number)
    props = (M = M, feh = MH(ts), vvcrit = ts.properties.vvcrit)
    nt = _generic_trackset_interp(ts, M)
    table = Table(NamedTuple{(:star_age, keys(nt)[2:end]...)}(tuple(exp10.(nt.logAge), values(nt)[2:end]...)))
    return MISTTrack(table, props)
end
mass(ts::MISTTrackSet) = ts.properties.masses
chemistry(::MISTTrackSet) = MISTChemistry()
MH(ts::MISTTrackSet) = ts.properties.feh
Z(ts::MISTTrackSet) = Z(chemistry(ts), MH(ts))
Y(ts::MISTTrackSet) = Y(chemistry(ts), Z(ts))
X(ts::MISTTrackSet) = 1 - Y(ts) - Z(ts)
post_rgb(t::MISTTrackSet) = true
Base.eltype(ts::MISTTrackSet) = typeof(ts.properties.feh)
function Base.show(io::IO, mime::MIME"text/plain", ts::MISTTrackSet)
    print(io, "MISTTrackSet with MH=$(MH(ts)), vvcrit=$(ts.properties.vvcrit), Z=$(Z(ts)), Y=$(Y(ts)), $(length(ts.AMRs)) EEPs and $(length(mass(ts))) initial stellar mass points.")
end

function isochrone(ts::MISTTrackSet, logAge::Number) # 1 ms
    eeps = Vector{Int}(undef, 0)
    track_extrema = extrema(mass(ts))
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

##########################################################################

"""
    MISTLibrary(vvcrit::Number=0)
`MISTLibrary` implements the [`AbstractTrackLibrary`](@ref StellarTracks.AbstractTrackLibrary)
interface for the MIST stellar evolution library. Instances can be constructed by providing a supported
`vvcrit` argument for the rotation parameter, which must be equal to either `0` (no rotation) or `0.4`.
We set `vvcrit=0` by default. If you construct an instance as `p = MISTLibrary(0.0)`, it is callable as
`p(mh::Number, M::Number)` which returns an [`InterpolatedTrack`](@ref StellarTracks.InterpolatedTrack)
that interpolates between tracks to a specific metallicity ([M/H]) and initial stellar mass (`M`).

This type also supports isochrone construction
(see [isochrone](@ref StellarTracks.isochrone(::StellarTracks.MIST.MISTLibrary, ::Number, ::Number))).

# Examples
```jldoctest
julia> p = MISTLibrary(0.0)
Structure of interpolants for the MIST library of stellar tracks with vvcrit=0.0. Valid range of metallicities is (-4.0, 0.5).

julia> isochrone(p, 10.05, -2) isa NamedTuple
true

julia> p(-2.05, 1.05)
InterpolatedTrack with M_ini=1.05, MH=-2.05, Z=0.0001327966689875739, Y=0.249199427865246, X=0.7506677754657665.
```
"""
struct MISTLibrary{A,B,C} <: AbstractTrackLibrary
    ts::A  # Vector of `TrackSet`s
    MH::B  # Vector of MH for each TrackSet
    vvcrit::C
end
chemistry(::MISTLibrary) = MISTChemistry()
MH(p::MISTLibrary) = p.MH # MH.(chemistry(tl), Z(tl))
Z(p::MISTLibrary) = Z.(chemistry(p), p.MH)
Y(p::MISTLibrary) = Y.(chemistry(p), Z(p))
X(p::MISTLibrary) = 1 .- Y(p) .- Z(p)
post_rgb(::MISTLibrary) = true
Base.eltype(p::MISTLibrary) = typeof(first(p.MH))
Base.Broadcast.broadcastable(p::MISTLibrary) = Ref(p)
function Base.show(io::IO, mime::MIME"text/plain", p::MISTLibrary)
    print(io, "Structure of interpolants for the MIST library of stellar tracks with vvcrit=$(p.vvcrit). Valid range of metallicities is $(extrema(MH(p))).")
end
function MISTLibrary(vvcrit::Number=0)
    # Make vector of tracksets
    ts = [MISTTrackSet(feh, vvcrit) for feh in feh_grid]
    return MISTLibrary(ts, feh_grid, vvcrit)
end

# Below is a stub for documentation,
# calls down to generic method in StellarTracks.jl main file
"""
    isochrone(p::MISTLibrary, logAge::Number, mh::Number)
Interpolates properties of the stellar tracks in the library at the requested logarithmic age (`logAge = log10(age [yr])`) and logarithmic metallicity [M/H] = `mh`. Returns a `NamedTuple` containing the properties listed below:
 - `eep`: Equivalent evolutionary points
 - `m_ini`: Initial stellar masses, in units of solar masses.
 - `logTe`: Base-10 logarithm of the effective temperature [K] of the stellar model.
 - `Mbol`: Bolometric luminosity of the stellar model.
 - `logg`: Surface gravity of the stellar model.
 - `log_surf_cell_z`: Base-10 logarithm of the surface metal mass fraction (Z).
"""
isochrone(p::MISTLibrary, logAge::Number, mh::Number)

export MISTTrack, MISTTrackSet, MISTLibrary, MISTChemistry   # Unique module exports
export mass, chemistry, X, Y, Z, MH, post_rgb, isochrone # Export generic API methods

end # module

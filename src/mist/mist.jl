"""StellarTracks.MIST provides access to the MIST stellar tracks."""
module MIST

# imports from parent module
using ..StellarTracks: AbstractTrack, AbstractTrackSet, AbstractTrackLibrary,
                       uniqueidx, Mbol, _generic_trackset_interp
import ..StellarTracks: X, Y, Z, MH, FeH, alphaFe, alpha_mass_fraction, chemistry, mass, post_rgb, isochrone, gridname # X_phot, Y_phot, Z_phot,

# Round floating-point values to 6 significant figures for display

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

# Reuse chemistry types from BolometricCorrections
using BolometricCorrections.MIST: MISTv1Chemistry, MISTv2Chemistry, unpack_txz

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
"""Initial stellar masses for the stellar tracks in the MIST v1.2 grid. These are uniform for all [Fe/H] and also for rotating and non-rotating grids."""
const mass_grid = track_type[0.1, 0.15, 0.2, 0.25, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.36, 1.38, 1.4, 1.42, 1.44, 1.46, 1.48, 1.5, 1.52, 1.54, 1.56, 1.58, 1.6, 1.62, 1.64, 1.66, 1.68, 1.7, 1.72, 1.74, 1.76, 1.78, 1.8, 1.82, 1.84, 1.86, 1.88, 1.9, 1.92, 1.94, 1.96, 1.98, 2.0, 2.02, 2.04, 2.06, 2.08, 2.1, 2.12, 2.14, 2.16, 2.18, 2.2, 2.22, 2.24, 2.26, 2.28, 2.3, 2.32, 2.34, 2.36, 2.38, 2.4, 2.42, 2.44, 2.46, 2.48, 2.5, 2.52, 2.54, 2.56, 2.58, 2.6, 2.62, 2.64, 2.66, 2.68, 2.7, 2.72, 2.74, 2.76, 2.78, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0]
"""Initial stellar masses for the stellar tracks in the MIST v2.5 grid. Note the gap between 9 M☉ and 40 M☉ for the [Fe/H] = -4.0, [α/Fe] = -0.2 grid; other metallicities may have different coverage."""
const mass_grid_v2 = track_type[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 40.0, 42.0, 44.0, 46.0, 48.0, 50.0, 52.0, 54.0, 56.0, 58.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0]
"""Available [Fe/H] values in the MIST v2.5 stellar track grid."""
const feh_grid_v2 = track_type[-4.0, -3.5, -3.0, -2.75, -2.5, -2.25, -2.0, -1.75, -1.5,
                               -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5]
"""Available [α/Fe] values in the MIST v2.5 stellar track grid."""
const afe_grid_v2 = track_type[-0.2, 0.0, 0.2, 0.4, 0.6]
"""Available vvcrit values in the MIST v2.5 stellar track grid."""
const vvcrit_grid_v2 = track_type[0.0, 0.4]

"""
    feh_grid_v2_for(afe) -> Vector
Return the subset of `feh_grid_v2` that is available for the given `[α/Fe]` value.
For `[α/Fe] = +0.6` the `[Fe/H] = +0.5` grid point is not available on the MIST server,
so it is excluded. All other `[α/Fe]` values have the full 17-point `feh_grid_v2`.
"""
function feh_grid_v2_for(afe::Number)
    afe_val = track_type(afe)
    if afe_val ≈ track_type(0.6)
        return filter(!≈(track_type(0.5)), feh_grid_v2)
    end
    return feh_grid_v2
end

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
    @argcheck !isnothing(mass_idx) ArgumentError("Invalid mass=$mass argument; available track masses are $mass_grid. For track interpolation, use `MISTv1Library`.")
    # mass = string(Int(mass * 100))
    mass = string(round(Int, mass_grid[mass_idx] * 100))
    mass = repeat("0", 5 - length(mass)) * mass # Pad to length 5
    return mass
end

"""
    _parse_mass_v2(mass::Number)
Given a MIST v2.5 initial stellar mass, validate and return its zero-padded filename
string (e.g. `1.0` → `"00100"`).
"""
function _parse_mass_v2(mass::Number)
    mass_idx = findfirst(≈(mass), mass_grid_v2)
    @argcheck !isnothing(mass_idx) ArgumentError("Invalid mass=$mass argument; available track masses are $mass_grid_v2. For track interpolation, use `MISTv2Library`.")
    mass_str = string(round(Int, mass_grid_v2[mass_idx] * 100))
    return repeat("0", 5 - length(mass_str)) * mass_str
end

"""
    _parse_afe_v2(afe::Number)
Validate `afe` against `afe_grid_v2` and return the canonical grid value.
"""
function _parse_afe_v2(afe::Number)
    afe_idx = findfirst(≈(afe), afe_grid_v2)
    @argcheck !isnothing(afe_idx) ArgumentError("Invalid afe=$afe argument; valid [α/Fe] values are $afe_grid_v2.")
    return afe_grid_v2[afe_idx]
end

##########################################################################

# Data download, organization, etc.
include("init.jl")

##########################################################################

"""
    MISTv1Track(mh::Number, mass::Number, vvcrit::Number=0)
`MISTv1Track` implements the [`AbstractTrack`](@ref StellarTracks.AbstractTrack)
interface for the MIST stellar evolution library.
```jldoctest
julia> track = StellarTracks.MIST.MISTv1Track(-2, 0.15, 0.0)
MISTv1Track with M_ini=0.15, MH=-2.0, vvcrit=0.0, Z=0.000148992, Y=0.249224, X=0.750627.

julia> track(7.0) # interpolate track at log10(age [yr]) = 7
(log_L = -1.5293719450743, log_Teff = 3.587337261741102, log_g = 4.447603584617846, log_surf_cell_z = -3.8450984758441953)
```
"""
struct MISTv1Track{A,B,C} <: AbstractTrack
    data::Table{A}
    itp::B
    properties::C
end
function MISTv1Track(data::Table, props)
    # Construct interpolator as a function of proper age
    # itp = interpolate(deduplicate_knots!(data.star_age; move_knots=true),
    #                   [SVector(values(i)[2:end]) for i in data], Gridded(Linear()))
    itp = CubicSpline([SVector(values(i)[2:end]) for i in data],
                      deduplicate_knots!(data.star_age; move_knots=true);
                      cache_parameters=true)
    return MISTv1Track(data, itp, props)
end
# Given feh, mass, vvcrit, load the data table and call to function above
function MISTv1Track(feh::Number, mass::Number, vvcrit::Number=0)
    props = (M = mass, feh = feh, vvcrit = vvcrit)
    feh_idx = findfirst(≈(feh), feh_grid) # Validate against feh_grid
    @argcheck !isnothing(feh_idx) ArgumentError("Provided `feh` argument $feh to `MISTv1Track` is invalid; available metallicities are $feh_grid. For metallicity interpolation, use `MISTv1Library`.")
    feh = string(feh_grid[feh_idx])
    # Validate vvcrit
    vvcrit = _parse_vvcrit(vvcrit)
    # Validate mass
    mass = _parse_mass(mass)
    # Load data file into table
    filename = @datadep_str(joinpath("MISTv1.2_vvcrit"*vvcrit, feh, mass * "M.track.jld2"))
    data = JLD2.load_object(filename)
    return MISTv1Track(data, props)
end
# Make Track callable with logAge to get properties as a NamedTuple
Base.keys(::MISTv1Track) = select_columns[2:end]
function (track::MISTv1Track)(logAge::Number)
    # result = track.itp(logAge)
    result = track.itp(exp10(logAge))
    return NamedTuple{keys(track)}(result)
end
Base.extrema(t::MISTv1Track) = log10.(extrema(t.itp.t))
gridname(::Type{<:MISTv1Track}) = "MISTv1"
mass(t::MISTv1Track) = t.properties.M
chemistry(::MISTv1Track) = MISTv1Chemistry()
MH(t::MISTv1Track) = t.properties.feh
FeH(t::MISTv1Track) = t.properties.feh
alphaFe(t::MISTv1Track) = zero(t.properties.feh)
# Whether there is post-rgb evolution or not is dependent on how many EEPs in track
post_rgb(t::MISTv1Track) = length(t.itp.u) > eep_idxs.RG_TIP
Base.eltype(t::MISTv1Track) = typeof(t.properties.feh)
function Base.show(io::IO, mime::MIME"text/plain", t::MISTv1Track)
    print(io, "MISTv1Track with M_ini=$(mass(t)), MH=$(round(MH(t); sigdigits=6)), vvcrit=$(t.properties.vvcrit), Z=$(round(Z(t); sigdigits=6)), Y=$(round(Y(t); sigdigits=6)), X=$(round(X(t); sigdigits=6)).")
end

##########################################################################

"""
    MISTv1TrackSet(mh::Number, vvcrit::Number=0)
`MISTv1TrackSet` implements the [`AbstractTrackSet`](@ref StellarTracks.AbstractTrackSet)
interface for the MIST stellar evolution library.
```jldoctest
julia> ts = StellarTracks.MIST.MISTv1TrackSet(0.0, 0.0)
MISTv1TrackSet with MH=0.0, vvcrit=0.0, Z=0.0142014, Y=0.270327, 1710 EEPs and 196 initial stellar mass points.

julia> ts(1.01) # Interpolate track at new initial mass
MISTv1Track with M_ini=1.01, MH=0.0, vvcrit=0.0, Z=0.0142014, Y=0.270327, X=0.715472.

julia> isochrone(ts, 10.0) isa NamedTuple # Interpolate isochrone at `log10(age [yr]) = 10`
true
```
"""
struct MISTv1TrackSet{A <: AbstractVector{<:Integer},
                    B <: AbstractVector{<:AbstractInterpolation},
                    C,
                    D} <: AbstractTrackSet
    eeps::A # EEP points; vector of integers
    AMRs::B # age-mass relations; vector of interpolators, same length as eeps
    interps::C
    properties::D
end
# Given a metallicity and rotation, load the correct models and then call below method
function MISTv1TrackSet(feh::Number, vvcrit::Number=0) # One table per stellar model
    # Validate feh
    feh_idx = findfirst(≈(feh), feh_grid) # Validate against feh_grid
    @argcheck !isnothing(feh_idx) ArgumentError("Provided `feh` argument $feh to `MISTv1TrackSet` is invalid; available metallicities are $feh_grid. For metallicity interpolation, use `MISTv1Library`.")
    feh = string(feh_grid[feh_idx])
    # Validate vvcrit
    vvcrit = _parse_vvcrit(vvcrit)
    # List stellar track files
    allfiles = [@datadep_str(joinpath("MISTv1.2_vvcrit"*vvcrit, feh, _parse_mass(mass) * "M.track.jld2")) for mass in mass_grid]
    # return [JLD2.load_object(file) for file in allfiles]
    # return MISTv1TrackSet([JLD2.load_object(file) for file in allfiles],
    #                     masses, parse(track_type, feh))
    data = vcat([begin
                     tmpdata = JLD2.load_object(allfiles[i])
                     Table(tmpdata, m_ini=fill(mass_grid[i], length(tmpdata)),
                           eep = 1:length(tmpdata))
                 end for i in eachindex(allfiles)]...)
    # return data
    return MISTv1TrackSet(data, parse(track_type, feh), parse(track_type, vvcrit))
end
function MISTv1TrackSet(data::Table, feh::Number, vvcrit::Number)
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
    return MISTv1TrackSet(eeps, amrs,
                        (log_L = logl, log_Teff = logte, log_g = logg, log_surf_cell_z = logsurfz),
                        (feh = feh, vvcrit = vvcrit, masses = unique(data.m_ini)))
end
function (ts::MISTv1TrackSet)(M::Number)
    props = (M = M, feh = MH(ts), vvcrit = ts.properties.vvcrit)
    nt = _generic_trackset_interp(ts, M)
    table = Table(NamedTuple{(:star_age, keys(nt)[2:end]...)}(tuple(exp10.(nt.logAge), values(nt)[2:end]...)))
    return MISTv1Track(table, props)
end
gridname(::Type{<:MISTv1TrackSet}) = "MISTv1"
mass(ts::MISTv1TrackSet) = ts.properties.masses
chemistry(::MISTv1TrackSet) = MISTv1Chemistry()
MH(ts::MISTv1TrackSet) = ts.properties.feh
FeH(ts::MISTv1TrackSet) = ts.properties.feh
alphaFe(ts::MISTv1TrackSet) = zero(ts.properties.feh)
post_rgb(t::MISTv1TrackSet) = true
Base.eltype(ts::MISTv1TrackSet) = typeof(ts.properties.feh)
function Base.show(io::IO, mime::MIME"text/plain", ts::MISTv1TrackSet)
    print(io, "MISTv1TrackSet with MH=$(round(MH(ts); sigdigits=6)), vvcrit=$(ts.properties.vvcrit), Z=$(round(Z(ts); sigdigits=6)), Y=$(round(Y(ts); sigdigits=6)), $(length(ts.AMRs)) EEPs and $(length(mass(ts))) initial stellar mass points.")
end

function isochrone(ts::MISTv1TrackSet, logAge::Number) # 1 ms
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
    MISTv1Library(vvcrit::Number=0)
`MISTv1Library` implements the [`AbstractTrackLibrary`](@ref StellarTracks.AbstractTrackLibrary)
interface for the MIST stellar evolution library. Instances can be constructed by providing a supported
`vvcrit` argument for the rotation parameter, which must be equal to either `0` (no rotation) or `0.4`.
We set `vvcrit=0` by default. If you construct an instance as `p = MISTv1Library(0.0)`, it is callable as
`p(mh::Number, M::Number)` which returns an [`InterpolatedTrack`](@ref StellarTracks.InterpolatedTrack)
that interpolates between tracks to a specific metallicity ([M/H]) and initial stellar mass (`M`).

This type also supports isochrone construction
(see [isochrone](@ref StellarTracks.isochrone(::StellarTracks.MIST.MISTv1Library, ::Number, ::Number))).

# Examples
```jldoctest
julia> p = MISTv1Library(0.0)
Structure of interpolants for the MIST library of stellar tracks with vvcrit=0.0. Valid range of metallicities is (-4.0, 0.5).

julia> isochrone(p, 10.05, -2) isa NamedTuple
true

julia> p(-2.05, 1.05)
InterpolatedTrack with M_ini=1.05, MH=-2.05, Z=0.000132797, Y=0.249199, X=0.750668.
```
"""
struct MISTv1Library{A,B,C} <: AbstractTrackLibrary
    ts::A  # Vector of `TrackSet`s
    MH::B  # Vector of MH for each TrackSet
    vvcrit::C
end
gridname(::Type{<:MISTv1Library}) = "MISTv1"
chemistry(::MISTv1Library) = MISTv1Chemistry()
MH(p::MISTv1Library) = p.MH # MH.(chemistry(tl), Z(tl))
FeH(p::MISTv1Library) = p.MH
alphaFe(p::MISTv1Library) = zero(eltype(p.MH))
post_rgb(::MISTv1Library) = true
Base.eltype(p::MISTv1Library) = typeof(first(p.MH))
Base.Broadcast.broadcastable(p::MISTv1Library) = Ref(p)
function Base.show(io::IO, mime::MIME"text/plain", p::MISTv1Library)
    print(io, "Structure of interpolants for the MIST library of stellar tracks with vvcrit=$(p.vvcrit). Valid range of metallicities is $(extrema(MH(p))).")
end
function MISTv1Library(vvcrit::Number=0)
    # Make vector of tracksets
    ts = [MISTv1TrackSet(feh, vvcrit) for feh in feh_grid]
    return MISTv1Library(ts, feh_grid, vvcrit)
end

# Below is a stub for documentation,
# calls down to generic method in StellarTracks.jl main file
"""
    isochrone(p::MISTv1Library, logAge::Number, mh::Number)
Interpolates properties of the stellar tracks in the library at the requested logarithmic age (`logAge = log10(age [yr])`) and logarithmic metallicity [M/H] = `mh`. Returns a `NamedTuple` containing the properties listed below:
 - `eep`: Equivalent evolutionary points
 - `m_ini`: Initial stellar masses, in units of solar masses.
 - `logTe`: Base-10 logarithm of the effective temperature [K] of the stellar model.
 - `Mbol`: Bolometric luminosity of the stellar model.
 - `logg`: Surface gravity of the stellar model.
 - `log_surf_cell_z`: Base-10 logarithm of the surface metal mass fraction (Z).
"""
isochrone(p::MISTv1Library, logAge::Number, mh::Number)

##########################################################################
# MIST v2.5 types
##########################################################################

"""
    MISTv2Track(feh, mass, vvcrit=0, afe=0)
`MISTv2Track` implements the [`AbstractTrack`](@ref StellarTracks.AbstractTrack)
interface for the MIST v2.5 stellar evolution library.
```jldoctest
julia> track = StellarTracks.MIST.MISTv2Track(-2, 0.15, 0.0, 0.0)
MISTv2Track with M_ini=0.15, MH=-2.0, afe=0.0, vvcrit=0.0, Z=0.000196117, Y=0.24926, X=0.750544.
```
"""
struct MISTv2Track{A,B,C} <: AbstractTrack
    data::Table{A}
    itp::B
    properties::C
end
function MISTv2Track(data::Table, props)
    itp = CubicSpline([SVector(values(i)[2:end]) for i in data],
                      deduplicate_knots!(data.star_age; move_knots=true);
                      cache_parameters=true)
    return MISTv2Track(data, itp, props)
end
function MISTv2Track(feh::Number, mass::Number, vvcrit::Number=0, afe::Number=0)
    props = (M = mass, feh = feh, vvcrit = vvcrit, afe = afe)
    valid_feh  = feh_grid_v2_for(afe)
    feh_idx = findfirst(≈(feh), valid_feh)
    @argcheck !isnothing(feh_idx) ArgumentError("Provided `feh` argument $feh to `MISTv2Track` is invalid for [α/Fe]=$afe; available metallicities are $valid_feh. For metallicity interpolation, use `MISTv2Library`.")
    feh_str    = string(valid_feh[feh_idx])
    vvcrit_str = _parse_vvcrit(vvcrit)
    afe_val    = _parse_afe_v2(afe)
    afe_str    = _afe_tag(afe_val)
    mass_str   = _parse_mass_v2(mass)
    filename = @datadep_str(joinpath("MISTv2.5_vvcrit$(vvcrit_str)_afe_$(afe_str)",
                                     feh_str, mass_str * "M.track.jld2"))
    data = JLD2.load_object(filename)
    return MISTv2Track(data, props)
end
Base.keys(::MISTv2Track) = select_columns[2:end]
function (track::MISTv2Track)(logAge::Number)
    result = track.itp(exp10(logAge))
    return NamedTuple{keys(track)}(result)
end
Base.extrema(t::MISTv2Track) = log10.(extrema(t.itp.t))
gridname(::Type{<:MISTv2Track}) = "MISTv2"
mass(t::MISTv2Track) = t.properties.M
chemistry(::MISTv2Track) = MISTv2Chemistry()
FeH(t::MISTv2Track) = t.properties.feh
alphaFe(t::MISTv2Track) = t.properties.afe
post_rgb(t::MISTv2Track) = length(t.itp.u) > eep_idxs.RG_TIP
Base.eltype(t::MISTv2Track) = typeof(t.properties.feh)
function Base.show(io::IO, mime::MIME"text/plain", t::MISTv2Track)
    print(io, "MISTv2Track with M_ini=$(mass(t)), MH=$(round(MH(t); sigdigits=6)), afe=$(t.properties.afe), vvcrit=$(t.properties.vvcrit), Z=$(round(Z(t); sigdigits=6)), Y=$(round(Y(t); sigdigits=6)), X=$(round(X(t); sigdigits=6)).")
end

##########################################################################

"""
    MISTv2TrackSet(feh, vvcrit=0, afe=0)
`MISTv2TrackSet` implements the [`AbstractTrackSet`](@ref StellarTracks.AbstractTrackSet)
interface for the MIST v2.5 stellar evolution library.
```jldoctest
julia> ts = StellarTracks.MIST.MISTv2TrackSet(0.0, 0.0, 0.0)
MISTv2TrackSet with MH=0.0, afe=0.0, vvcrit=0.0, Z=0.0185, Y=0.2735, 1721 EEPs and 155 initial stellar mass points.
```
"""
struct MISTv2TrackSet{A <: AbstractVector{<:Integer},
                      B <: AbstractVector{<:AbstractInterpolation},
                      C,
                      D} <: AbstractTrackSet
    eeps::A
    AMRs::B
    interps::C
    properties::D
end
function MISTv2TrackSet(feh::Number, vvcrit::Number=0, afe::Number=0)
    valid_feh  = feh_grid_v2_for(afe)
    feh_idx = findfirst(≈(feh), valid_feh)
    @argcheck !isnothing(feh_idx) ArgumentError("Provided `feh` argument $feh to `MISTv2TrackSet` is invalid for [α/Fe]=$afe; available metallicities are $valid_feh. For metallicity interpolation, use `MISTv2Library`.")
    feh_str    = string(valid_feh[feh_idx])
    vvcrit_str = _parse_vvcrit(vvcrit)
    afe_val    = _parse_afe_v2(afe)
    afe_str    = _afe_tag(afe_val)
    datadep_dir = @datadep_str("MISTv2.5_vvcrit$(vvcrit_str)_afe_$(afe_str)")
    allfiles_masses = [(path = joinpath(datadep_dir, feh_str, _parse_mass_v2(m) * "M.track.jld2"),
                        mass = m)
                       for m in mass_grid_v2]
    # Some (feh, vvcrit, afe) combinations are missing individual mass tracks in the MIST v2.5 grid.
    # Skip files that are not present on disk rather than erroring.
    allfiles_masses = filter(fm -> isfile(fm.path), allfiles_masses)
    data = vcat([begin
                     tmpdata = JLD2.load_object(fm.path)
                     Table(tmpdata, m_ini=fill(fm.mass, length(tmpdata)),
                           eep = 1:length(tmpdata))
                 end for fm in allfiles_masses]...)
    return MISTv2TrackSet(data, valid_feh[feh_idx],
                          parse(track_type, vvcrit_str), afe_val)
end
function MISTv2TrackSet(data::Table, feh::Number, vvcrit::Number, afe::Number)
    eeps = sort(unique(data.eep))
    itp_type = CubicHermiteSpline{Vector{track_type},
                                  Vector{track_type},
                                  Vector{track_type},
                                  Vector{track_type},
                                  Vector{track_type},
                                  track_type}
    amrs    = Vector{itp_type}(undef, length(eeps))
    logte   = Vector{itp_type}(undef, length(eeps))
    logl    = similar(logte)
    logg    = similar(logte)
    logsurfz = similar(logte)

    Threads.@threads for i in eachindex(eeps)
        eep = eeps[i]
        tmpdata = data[findall(Base.Fix1(==, eep), data.eep)]
        idxs = sortperm(tmpdata.star_age)
        tmpdata = tmpdata[idxs]
        tmpdata = tmpdata[uniqueidx(tmpdata.star_age)]
        _, p1 = findmax(tmpdata.m_ini)
        idxs = p1:lastindex(tmpdata)
        goodidxs = vcat(true, diff(tmpdata.m_ini[idxs]) .< 0)
        tmpdata = tmpdata[idxs[goodidxs]]
        amrs[i] = PCHIPInterpolation(tmpdata.m_ini, log10.(tmpdata.star_age))
        idxs = sortperm(tmpdata.m_ini)
        tmpdata = tmpdata[idxs]
        logte[i]    = PCHIPInterpolation(tmpdata.log_Teff, tmpdata.m_ini)
        logl[i]     = PCHIPInterpolation(tmpdata.log_L, tmpdata.m_ini)
        logg[i]     = PCHIPInterpolation(tmpdata.log_g, tmpdata.m_ini)
        logsurfz[i] = PCHIPInterpolation(tmpdata.log_surf_cell_z, tmpdata.m_ini)
    end
    return MISTv2TrackSet(eeps, amrs,
                         (log_L = logl, log_Teff = logte, log_g = logg, log_surf_cell_z = logsurfz),
                         (feh = feh, vvcrit = vvcrit, afe = afe, masses = unique(data.m_ini)))
end
function (ts::MISTv2TrackSet)(M::Number)
    props = (M = M, feh = FeH(ts), vvcrit = ts.properties.vvcrit, afe = ts.properties.afe)
    nt = _generic_trackset_interp(ts, M)
    table = Table(NamedTuple{(:star_age, keys(nt)[2:end]...)}(
                  tuple(exp10.(nt.logAge), values(nt)[2:end]...)))
    return MISTv2Track(table, props)
end
gridname(::Type{<:MISTv2TrackSet}) = "MISTv2"
mass(ts::MISTv2TrackSet) = ts.properties.masses
chemistry(::MISTv2TrackSet) = MISTv2Chemistry()
FeH(ts::MISTv2TrackSet) = ts.properties.feh
alphaFe(ts::MISTv2TrackSet) = ts.properties.afe
post_rgb(::MISTv2TrackSet) = true
Base.eltype(ts::MISTv2TrackSet) = typeof(ts.properties.feh)
function Base.show(io::IO, mime::MIME"text/plain", ts::MISTv2TrackSet)
    print(io, "MISTv2TrackSet with MH=$(round(MH(ts); sigdigits=6)), afe=$(ts.properties.afe), vvcrit=$(ts.properties.vvcrit), Z=$(round(Z(ts); sigdigits=6)), Y=$(round(Y(ts); sigdigits=6)), $(length(ts.AMRs)) EEPs and $(length(mass(ts))) initial stellar mass points.")
end

function isochrone(ts::MISTv2TrackSet, logAge::Number)
    eeps = Vector{Int}(undef, 0)
    track_extrema = extrema(mass(ts))
    interp_masses = Vector{eltype(ts)}(undef, 0)
    logte   = similar(interp_masses)
    logl    = similar(interp_masses)
    logg    = similar(interp_masses)
    logsurfz = similar(interp_masses)
    for (i, amr) in enumerate(ts.AMRs)
        drange = extrema(amr)
        if logAge >= first(drange) && logAge <= last(drange)
            imass = amr(logAge)
            if imass >= first(track_extrema) && imass <= last(track_extrema)
                logli = ts.interps.log_L[i](imass)
                if length(logl) > 0 && i < eep_idxs[4] && logli < last(logl)
                    li = lastindex(logl)
                    deleteat!(eeps, li); deleteat!(interp_masses, li)
                    deleteat!(logte, li); deleteat!(logl, li)
                    deleteat!(logg, li); deleteat!(logsurfz, li)
                    continue
                end
                push!(eeps, i); push!(interp_masses, imass)
                push!(logl, logli)
                push!(logte,    ts.interps.log_Teff[i](imass))
                push!(logg,     ts.interps.log_g[i](imass))
                push!(logsurfz, ts.interps.log_surf_cell_z[i](imass))
            end
        end
    end
    return (eep = eeps, m_ini = interp_masses, logTe = logte, Mbol = Mbol.(logl, 4.74),
            logg = logg, logL = logl, log_surf_cell_z = logsurfz)
end

##########################################################################

"""
    MISTv2Library(vvcrit=0.0, afe=0.0)
`MISTv2Library` implements the [`AbstractTrackLibrary`](@ref StellarTracks.AbstractTrackLibrary)
interface for the MIST v2.5 stellar evolution library.

The library interpolates in [M/H] (metallicity). The `vvcrit` and `afe` parameters
must be exact values from the MIST v2.5 grid (`vvcrit_grid_v2` and `afe_grid_v2`);
they select which DataDep is loaded but are not interpolated over.

Calling the library as `lib(mh, mass)` returns a
[`MISTv2Track`](@ref) or [`InterpolatedTrack`](@ref StellarTracks.InterpolatedTrack)
interpolated to `mh` and `mass`.

This type also supports isochrone construction
(see [isochrone](@ref isochrone(::MISTv2Library, ::Number, ::Number))).

# Examples
```jldoctest
julia> p = MISTv2Library(0.0, 0.0)
MISTv2Library with vvcrit=0.0, afe=0.0. Valid [M/H] range: (-4.0, 0.5).

julia> isochrone(p, 10.05, -2) isa NamedTuple
true

julia> p(-2.05, 1.05)
InterpolatedTrack with M_ini=1.05, MH=-2.05, Z=0.000174801, Y=0.249231, X=0.750594.
```
"""
struct MISTv2Library{A, B, C, D} <: AbstractTrackLibrary
    ts::A  # Vector{MISTv2TrackSet}
    feh::B
    vvcrit::C
    afe::D
end
gridname(::Type{<:MISTv2Library}) = "MISTv2"
chemistry(::MISTv2Library) = MISTv2Chemistry()
FeH(p::MISTv2Library) = p.feh
alphaFe(p::MISTv2Library) = p.afe
post_rgb(::MISTv2Library) = true
Base.eltype(p::MISTv2Library) = eltype(p.feh)
Base.Broadcast.broadcastable(p::MISTv2Library) = Ref(p)
function Base.show(io::IO, mime::MIME"text/plain", p::MISTv2Library)
    print(io, "MISTv2Library with vvcrit=$(p.vvcrit), afe=$(p.afe). ",
          "Valid [M/H] range: $(extrema(MH(p))).")
end

function MISTv2Library(vvcrit::Number=0.0, afe::Number=0.0)
    vvcrit_val = track_type(vvcrit)
    afe_val    = track_type(afe)
    @argcheck any(≈(vvcrit_val), vvcrit_grid_v2) ArgumentError("Invalid vvcrit=$vvcrit; valid options are $vvcrit_grid_v2.")
    @argcheck any(≈(afe_val),    afe_grid_v2)    ArgumentError("Invalid afe=$afe; valid options are $afe_grid_v2.")
    valid_feh = feh_grid_v2_for(afe_val)
    ts = [MISTv2TrackSet(feh, vvcrit_val, afe_val) for feh in valid_feh]
    return MISTv2Library(ts, valid_feh, vvcrit_val, afe_val)
end

# Below is a stub for documentation,
# calls down to generic method in StellarTracks.jl main file
"""
    isochrone(p::MISTv2Library, logAge::Number, mh::Number)
Interpolates properties of the stellar tracks in the library at the requested logarithmic age (`logAge = log10(age [yr])`) and logarithmic metallicity [M/H] = `mh`. Returns a `NamedTuple` containing the properties listed below:
 - `eep`: Equivalent evolutionary points
 - `m_ini`: Initial stellar masses, in units of solar masses.
 - `logTe`: Base-10 logarithm of the effective temperature [K] of the stellar model.
 - `Mbol`: Bolometric luminosity of the stellar model.
 - `logg`: Surface gravity of the stellar model.
 - `log_surf_cell_z`: Base-10 logarithm of the surface metal mass fraction (Z).
"""
isochrone(p::MISTv2Library, logAge::Number, mh::Number)

export MISTv1Track, MISTv1TrackSet, MISTv1Library, MISTv1Chemistry   # Unique module exports
export MISTv2Track, MISTv2TrackSet, MISTv2Library, MISTv2Chemistry
export mass, chemistry, X, Y, Z, MH, FeH, alphaFe, alpha_mass_fraction, post_rgb, isochrone, gridname # Export generic API methods

# Deprecated aliases — forward to v1 types
@deprecate MISTTrack(args...) MISTv1Track(args...)
@deprecate MISTTrackSet(args...) MISTv1TrackSet(args...)
@deprecate MISTLibrary(args...) MISTv1Library(args...)
export MISTTrack, MISTTrackSet, MISTLibrary  # keep exported so old code can use them

end # module

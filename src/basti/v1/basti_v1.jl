"""StellarTracks.BASTIv1 provides access to the older BaSTI stellar tracks (circa 2013)."""
module BASTIv1

# imports from parent module
using ..StellarTracks: AbstractChemicalMixture, AbstractTrack, AbstractTrackSet,
                       AbstractTrackLibrary, uniqueidx, Mbol, _generic_trackset_interp
import ..StellarTracks: X, Y, Z, X_phot, Y_phot, Z_phot, MH, chemistry, mass, post_rgb, isochrone

# Imports for data reading / processing
# import CSV
using DataDeps: register, DataDep, @datadep_str
using DataInterpolations: AbstractInterpolation, CubicSpline, CubicHermiteSpline, PCHIPInterpolation
# using DelimitedFiles: readdlm
# using Glob: glob
import JLD2 # for saving files in binary format
using Printf: @sprintf
using TypedTables: Table

# Imports for core module code
using ArgCheck: @argcheck
using StaticArrays: SVector

"""Number of secondary EEP points per primary EEP point. Pietrinferni 2004 Table 2, called key points (KPs). """
const eep_lengths = (MS_BEG = 200,    # beginning of MS
                     # intermediate age MS; min(Teff) for high-mass (HM); X_c = 0.3 for low-mass (LM)
                     IAMS = 60,       
                     # MS_KINK = 61,    # MS kink point; maximum Teff along MS
                     MS_TO = 60,      # terminal MS; max(L) for HM or X_c = 0 for LM
                     SGB = 70,        # MSTO to end of subgiant branch; min(L) for HM
                     RG_START = 800,  # RGB start to RGB tip
                     HE_BEG = 10,     # beginning core He fusion
                     HE_2 = 150,      # core He fusion Y_c = 0.55
                     HE_3 = 100,      # core He fusion Y_c = 0.50
                     HE_4 = 100,      # core He fusion Y_c = 0.40
                     HE_5 = 80,       # core He fusion Y_c = 0.20
                     HE_6 = 80,       # core He fusion Y_c = 0.10
                     END_CHEB = 140,  # end of core He burning; core He fusion Y_c = 0.0
                     TPAGB_BEG = 150) # L_{CNO} > L_{3α} during AGB phase
# Indices where the different phases begin
"""1-based indices giving the EEP point at which each EEP phase begins."""
const eep_idxs = NamedTuple{keys(eep_lengths)}((1, (cumsum(values(eep_lengths)[begin:end-1]) .+ 1)...))
"""Data type to parse the BaSTIv1 tracks as."""
const track_type = Float32
"""Valid metal mass fractions (Z) for BaSTIv1."""
const zgrid = track_type[0.00001, 0.0001, 0.0003, 0.0006, 0.001, 0.002, 0.004, 0.008, 0.01, 0.0198, 0.003, 0.04, 0.05]
"""Valid helium mass fractions (Z) for BaSTIv1."""
const ygrid = track_type[0.245, 0.245, 0.245, 0.246, 0.246, 0.248, 0.251, 0.256, 0.259, 0.2734, 0.288, 0.303, 0.316]
"""Initial stellar masses for the stellar tracks in the BaSTIv1 grid. These are uniform for all metallicities and also for alpha-enhanced grids."""
const massgrid = track_type[0.5, 0.55, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.45, 2.5, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
"""α-element enrichment parameters [α/Fe] available for the BaSTIv1 grid."""
const αFegrid = track_type[0.0, 0.4]

"""
    _parse_zval(zval::Number)
Convert a metal mass fraction `zval` into a properly formatted string for use in loading the JLD2 file.
"""
_parse_zval(zval::Number) = @sprintf("%.5f", zval)
"""
    _parse_zval(_parse_α_fe(α_fe::Number))
Convert an alpha enrichment factor `α_fe` into a properly formatted string for use in loading the JLD2 file.
"""
_parse_α_fe(α_fe::Number) = @sprintf("%.1f", α_fe)

##########################################################################
# BaSTIv1 chemistry
"""
    BaSTIv1Chemistry()
Returns a singleton struct representing the BaSTIv1 chemical mixture model.
These older BaSTI models, presented in [Pietrinferni2004,Pietrinferni2006,Pietrinferni2013](@citet), include both solar-scaled chemical compositions and α-enhanced compositions with
[α/Fe] ≈ 0.4. The solar protostellar chemical mixture for these models was calibrated to
reproduce solar photospheric observations via a forward modeling approach
(see section 4 of [Pietrinferni2004](@citet)).
The distribution of heavy metals is taken from [Grevesse1993](@citet), which the authors
state is minimally different from [Grevesse1998](@citet).

```jldoctest
julia> using StellarTracks.BASTIv1: BaSTIv1Chemistry, X, Y, Z, X_phot, Y_phot, Z_phot, MH;

julia> chem = BaSTIv1Chemistry();

julia> X(chem) + Y(chem) + Z(chem) ≈ 1 # solar protostellar values
true

julia> X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1 # solar photospheric values
true

julia> MH(chem, Z(chem) * 0.1) ≈ -1.025908305527913
true

julia> Z(chem, -1.025908305527913) ≈ Z(chem) * 0.1
true
```
"""
struct BaSTIv1Chemistry <: AbstractChemicalMixture end
X(::BaSTIv1Chemistry) = 0.7068 # Protostellar abundance of calibrated solar model
X_phot(::BaSTIv1Chemistry) = 0.73709 # Calculated from (Z/X) = 0.0244
Y(::BaSTIv1Chemistry) = 0.2734 # Protostellar abundance of calibrated solar model
Y_phot(::BaSTIv1Chemistry) = 0.244 # Assumed present-day photospheric abundance
Y_p(::BaSTIv1Chemistry) = 0.245  # Second paragraph, section 5.1, Pietrinferni2004
Z(::BaSTIv1Chemistry) = 0.0198 # Protostellar abundance of calibrated solar model
Z_phot(::BaSTIv1Chemistry) = 0.01891 # Calculated from (Z/X) = 0.0244

Y(mix::BaSTIv1Chemistry, Zval) = Y_p(mix) + 1.4 * Zval # Second paragraph, section 5.1, Pietrinferni2004
# X generic
MH(mix::BaSTIv1Chemistry, Zval) = log10(Zval / X(mix, Zval)) - log10(Z(mix) / X(mix))
function Z(mix::BaSTIv1Chemistry, MHval)
    # Derivation in parsec code
    zoverx = exp10(MHval + log10(Z(mix) / X(mix)))
    γ = 1.4
    return (1 - Y_p(mix)) * zoverx / (1 + (1 + γ) * zoverx)
end

##########################################################################

# Data download, organization, etc.
include("init.jl")

##########################################################################

# """
#     MISTTrack(mh::Number, mass::Number, vvcrit::Number=0)
# `MISTTrack` implements the [`AbstractTrack`](@ref StellarTracks.AbstractTrack)
# interface for the MIST stellar evolution library.
# ```jldoctest
# julia> track = StellarTracks.MIST.MISTTrack(-2, 0.15, 0.0)
# MISTTrack with M_ini=0.15, MH=-2, vvcrit=0.0, Z=0.00014899227095992976, Y=0.24922374966753474, X=0.7506272580615054.

# julia> track(7.0) # interpolate track at log10(age [yr]) = 7
# (log_L = -1.5293719450743, log_Teff = 3.587337261741102, log_g = 4.447603584617846, log_surf_cell_z = -3.8450984758441953)
# ```
# """
# struct MISTTrack{A,B,C} <: AbstractTrack
#     data::Table{A}
#     itp::B
#     properties::C
# end
# function MISTTrack(data::Table, props)
#     # Construct interpolator as a function of proper age
#     # itp = interpolate(deduplicate_knots!(data.star_age; move_knots=true),
#     #                   [SVector(values(i)[2:end]) for i in data], Gridded(Linear()))
#     itp = CubicSpline([SVector(values(i)[2:end]) for i in data],
#                       deduplicate_knots!(data.star_age; move_knots=true))
#     return MISTTrack(data, itp, props)
# end
# # Given feh, mass, vvcrit, load the data table and call to function above
# function MISTTrack(feh::Number, mass::Number, vvcrit::Number=0)
#     props = (M = mass, feh = feh, vvcrit = vvcrit)
#     # Validate feh
#     @argcheck feh in feh_grid
#     feh = string(feh_grid[searchsortedfirst(feh_grid, feh)])
#     # Validate vvcrit
#     vvcrit = _parse_vvcrit(vvcrit)
#     dd_path = @datadep_str("MISTv1.2_vvcrit"*vvcrit)
#     # Validate mass
#     allfiles = readdir(joinpath(dd_path, feh); join=true)
#     masses = mist_mass.(allfiles)
#     @argcheck mass in masses ArgumentError("Invalid mass=$mass argument; available track masses for [Fe/H]=$feh and vvcrit=$vvcrit are $masses.")
#     # mass = string(Int(mass * 100))
#     mass = string(round(Int, mass * 100))
#     mass = repeat("0", 5 - length(mass)) * mass # Pad to length 5
#     # Load data file into table
#     filename = joinpath(dd_path, feh, mass * "M.track.jld2")
#     data = JLD2.load_object(filename)
#     return MISTTrack(data, props)
# end
# # Make Track callable with logAge to get properties as a NamedTuple
# function (track::MISTTrack)(logAge::Number)
#     # result = track.itp(logAge)
#     result = track.itp(exp10(logAge))
#     return NamedTuple{Tuple(select_columns)[2:end]}(result)
# end
# Base.extrema(t::MISTTrack) = log10.(extrema(t.itp.t))
# mass(t::MISTTrack) = t.properties.M
# chemistry(::MISTTrack) = MISTChemistry()
# MH(t::MISTTrack) = t.properties.feh # MH(chemistry(t), Z(t))
# Z(t::MISTTrack) = Z(chemistry(t), MH(t)) # t.properties.Z
# Y(t::MISTTrack) = Y(chemistry(t), Z(t))  # t.properties.Y
# X(t::MISTTrack) = 1 - Y(t) - Z(t)
# # Whether there is post-rgb evolution or not is dependent on how many EEPs in track
# post_rgb(t::MISTTrack) = length(t.itp.u) > eep_idxs.RG_TIP
# Base.eltype(t::MISTTrack) = typeof(t.properties.feh)
# function Base.show(io::IO, mime::MIME"text/plain", t::MISTTrack)
#     print(io, "MISTTrack with M_ini=$(mass(t)), MH=$(MH(t)), vvcrit=$(t.properties.vvcrit), Z=$(Z(t)), Y=$(Y(t)), X=$(X(t)).")
# end

##########################################################################

"""
    BaSTIv1TrackSet(zval::Number, α_fe::Number=0, canonical::Bool=false)
`BaSTIv1TrackSet` implements the [`AbstractTrackSet`](@ref StellarTracks.AbstractTrackSet)
interface for the older BaSTI stellar evolution library.
```jldoctest
julia> ts = StellarTracks.BASTIv1.BaSTIv1TrackSet(1e-3, 0.0, false)
Non-canonical, BaSTIv1TrackSet with MH=-1.3239328427169887, [α/Fe]=0.0, Z=0.001, Y=0.24640000006649643, 1999 EEPs and 25 initial stellar mass points.

julia> ts(1.01) # Interpolate track at new initial mass

julia> isochrone(ts, 10.0) isa NamedTuple # Interpolate isochrone at `log10(age [yr]) = 10`
true
```
"""
struct BaSTIv1TrackSet{A <: AbstractVector{<:Integer},
                       B <: AbstractVector{<:AbstractInterpolation},
                       C,
                       D} <: AbstractTrackSet
    eeps::A # EEP points; vector of integers
    AMRs::B # age-mass relations; vector of interpolators, same length as eeps
    interps::C
    properties::D
end
# Given a metallicity, α-enhancement, and canonical status, load the correct models and then call below method
function BaSTIv1TrackSet(zval::Number, α_fe::Number=0, canonical::Bool=false)
    # Validate Z, [α/Fe]
    @argcheck any(Base.Fix1(isapprox, zval), zgrid)
    @argcheck any(Base.Fix1(isapprox, α_fe), αFegrid)
    # Convert to strings
    zval = _parse_zval(zgrid[searchsortedfirst(zgrid, zval)])
    α_fe = _parse_α_fe(αFegrid[searchsortedfirst(αFegrid, α_fe)])
    # Validate vvcrit
    dd_path = @datadep_str("BaSTIv1")
    group = ifelse(canonical, "canonical", "noncanonical") * "/" * α_fe * "/" * zval
    # println(group)
    bfile = joinpath(dd_path, "basti2013.jld2")
    data = JLD2.load(bfile, group)
    # eeps = reduce(vcat, [begin
    #                          tmpdata = data[data.m_ini .== m]
    #                          1:length(tmpdata)
    #                      end for m in unique(data.m_ini)])
    # data = Table(eeps = eeps, data)
    # return data
    return BaSTIv1TrackSet(data, parse(track_type, zval), parse(track_type, α_fe), canonical)
end
function BaSTIv1TrackSet(data::Table, zval::Number, α_fe::Number, canonical::Bool)
    eeps = sort(unique(data.eep))[2:end] # First EEP has same age for all masses, so skip
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

    for i in eachindex(eeps)
        eep = eeps[i]
        tmpdata = data[findall(Base.Fix1(==, eep), data.eep)] # Performance optimization
        # Sort by age, which will be the independent variable in the AMR interpolation
        idxs = sortperm(tmpdata.logAge)
        tmpdata = tmpdata[idxs]
        # Now keep only entries that define a unique mapping between age and m_ini
        tmpdata = tmpdata[uniqueidx(tmpdata.logAge)]

        # We'll start the interpolation at the most massive model and enforce a monotonic
        # decrease in stellar initial mass with increasing age at fixed EEP
        _, p1 = findmax(tmpdata.m_ini)
        idxs = p1:lastindex(tmpdata)
        goodidxs = diff(tmpdata.m_ini[idxs]) .< 0
        goodidxs = vcat(true, goodidxs) # add true for first element as well
        tmpdata = tmpdata[idxs[goodidxs]]

        # PCHIPInterpolation is a type of CubicHermiteSpline
        amrs[i] = PCHIPInterpolation(tmpdata.m_ini, tmpdata.logAge)
        # Sort by initial stellar mass for defining the other interpolations
        idxs = sortperm(tmpdata.m_ini)
        tmpdata = tmpdata[idxs]
        # Now interpolate mass against logte, mbol, logg, c_o
        logte[i] = PCHIPInterpolation(log10.(tmpdata.Teff), tmpdata.m_ini)
        logl[i] = PCHIPInterpolation(tmpdata.logL, tmpdata.m_ini)
        logg[i] = PCHIPInterpolation(tmpdata.logg, tmpdata.m_ini)
    end
    return BaSTIv1TrackSet(eeps, amrs,
                           (log_L = logl, log_Teff = logte, log_g = logg),
                           (Z = zval, α_fe = α_fe, canonical = canonical,
                           masses = unique(data.m_ini)))
end
# function (ts::BaSTIv1TrackSet)(M::Number)
#     props = (M = M, feh = MH(ts), vvcrit = ts.properties.vvcrit)
#     nt = _generic_trackset_interp(ts, M)
#     table = Table(NamedTuple{(:star_age, keys(nt)[2:end]...)}(tuple(exp10.(nt.logAge), values(nt)[2:end]...)))
#     return MISTTrack(table, props)
# end
mass(ts::BaSTIv1TrackSet) = ts.properties.masses
chemistry(::BaSTIv1TrackSet) = BaSTIv1Chemistry()
Z(ts::BaSTIv1TrackSet) = ts.properties.Z
MH(ts::BaSTIv1TrackSet) = MH(chemistry(ts), Z(ts))
Y(ts::BaSTIv1TrackSet) = Y(chemistry(ts), Z(ts))
X(ts::BaSTIv1TrackSet) = 1 - Y(ts) - Z(ts)
post_rgb(t::BaSTIv1TrackSet) = true
Base.eltype(ts::BaSTIv1TrackSet) = typeof(ts.properties.feh)
function Base.show(io::IO, mime::MIME"text/plain", ts::BaSTIv1TrackSet)
    print(io, """$(ifelse(ts.properties.canonical, "Canonical", "Non-canonical")) BaSTIv1TrackSet with MH=$(MH(ts)), [α/Fe]=$(ts.properties.α_fe), Z=$(Z(ts)), Y=$(Y(ts)), $(length(ts.AMRs)) EEPs and $(length(ts.properties.masses)) initial stellar mass points.""")
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

##########################################################################

"""
    MISTLibrary(vvcrit::Number=0)
`MISTLibrary` implements the [`AbstractTrackLibrary`](@ref StellarTracks.AbstractTrackLibrary)
interface for the MIST stellar evolution library. Instances can be constructed by providing a supported
`vvcrit` argument for the rotation parameter, which must be equal to either `0` (no rotation) or `0.4`.
We set `vvcrit=0` by default. If you construct an instance as `p = MISTLibrary(0.0)`, it is callable as
 - `p(mh::Number)` to interpolate the full library to a new metallicity
   (returning a [`MISTTrackSet`](@ref)), or
 - `p(mh::Number, M::Number)` to interpolate the tracks to a specific metallicity
   and initial stellar mass (returning a [`MISTTrack`](@ref)).

This type also supports isochrone construction
(see [isochrone](@ref StellarTracks.isochrone(::StellarTracks.MIST.MISTLibrary, ::Number, ::Number))).

# Examples
```jldoctest
julia> p = MISTLibrary(0.0)
Structure of interpolants for the MIST library of stellar tracks with vvcrit=0.0. Valid range of metallicities is (-4.0, 0.5).

julia> isochrone(p, 10.05, -2) isa NamedTuple
true
```
"""
struct MISTLibrary{A,B,C} <: AbstractTrackLibrary
    ts::A  # Vector of `TrackSet`s
    MH::B  # Vector of MH for each TrackSet
    vvcrit::C
end
# Interpolation to get a TrackSet with metallicity MH
function (ts::MISTLibrary)(mh::Number)
    error("Not yet implemented.")
end
# Interpolation to get a Track with mass M and metallicity MH
function (ts::MISTLibrary)(mh::Number, M::Number)
    error("Not yet implemented.")
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

# export PARSECLibrary, PARSECChemistry, MH_canon, Z_canon # Unique module exports
export MISTTrack, MISTTrackSet, MISTLibrary, MISTChemistry   # Unique module exports
export mass, chemistry, X, Y, Z, MH, post_rgb, isochrone # Export generic API methods

end # module

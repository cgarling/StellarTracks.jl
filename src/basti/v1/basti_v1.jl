"""StellarTracks.BASTIv1 provides access to the older BaSTI stellar tracks (circa 2013)."""
module BaSTIv1

# imports from parent module
using ..StellarTracks: AbstractChemicalMixture, AbstractTrack, AbstractTrackSet,
                       AbstractTrackLibrary, uniqueidx, Mbol, _generic_trackset_interp,
                       radius, surface_gravity
import ..StellarTracks: X, Y, Z, X_phot, Y_phot, Z_phot, MH, chemistry, mass, post_rgb, isochrone

# Imports for data reading / processing
using DataDeps: register, DataDep, @datadep_str
using DataInterpolations: AbstractInterpolation, CubicSpline, CubicHermiteSpline, PCHIPInterpolation
using Interpolations: deduplicate_knots!
import JLD2 # for saving files in binary format
using Printf: @sprintf
using TypedTables: Table, columnnames

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
const zgrid = track_type[0.00001, 0.0001, 0.0003, 0.0006, 0.001, 0.002, 0.004, 0.008, 0.01, 0.0198, 0.03, 0.04, 0.05]
"""Valid helium mass fractions (Z) for BaSTIv1."""
const ygrid = track_type[0.245, 0.245, 0.245, 0.246, 0.246, 0.248, 0.251, 0.256, 0.259, 0.2734, 0.288, 0.303, 0.316]
"""Initial stellar masses for the stellar tracks in the BaSTIv1 grid. These are uniform for all metallicities and also for alpha-enhanced grids."""
const massgrid = track_type[0.5, 0.55, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.45, 2.5, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
"""α-element enrichment parameters [α/Fe] available for the BaSTIv1 grid."""
const αFegrid = track_type[0.0, 0.4]
"""Values of the Reimers mass loss parameter available for the BaSTIv1 grid."""
const ηgrid = track_type[0.2, 0.4]

function _validate_params(zval::Number, α_fe::Number, canonical::Bool, agb::Bool, η::Number)
    if ~any(≈(zval), zgrid)
        throw(ArgumentError("Input Z $zval invalid; available Z values are $(zgrid)."))
    end
    if (1e-5 ≈ zval || 0.05 ≈ zval) && (~(η ≈ 0.4) || agb)
        throw(ArgumentError("""Z value $zval is only available for η=0.4 models without AGB extension; you requested η=$η $(ifelse(agb, "with AGB extension", "without AGB extension"))."""))
    end
    if ~any(≈(α_fe), αFegrid)
        throw(ArgumentError("Supported [α/Fe] values for the BaSTIv1 grid are 0.0, 0.4; you requested $α_fe."))
    end
    if ~any(≈(η), ηgrid)
        throw(ArgumentError("Supported values of the Reimers mass loss parameter η are 0.2, 0.4; you requested $η."))
    end
    return true
end

"""
    _parse_zval(zval::Number)
Convert a metal mass fraction `zval` into a properly formatted string for use in loading the JLD2 file.
"""
_parse_zval(zval::Number) = @sprintf("%1.5f", zval)
"""
    _parse_α_fe(α_fe::Number)
Convert an alpha enrichment factor `α_fe` into a properly formatted string for use in loading the JLD2 file.
"""
_parse_α_fe(α_fe::Number) = @sprintf("%1.1f", α_fe)
"""
    _parse_η(η::Number)
Convert an alpha enrichment factor `α_fe` into a properly formatted string for use in loading the JLD2 file.
"""
_parse_η(η::Number) = @sprintf("%1.1f", η)

##########################################################################
# BaSTIv1 chemistry
"""
    BaSTIv1Chemistry()
Returns a singleton struct representing the BaSTIv1 chemical mixture model.
These older BaSTI models, presented in [Pietrinferni2004,Pietrinferni2006,Pietrinferni2013](@citet),
include both solar-scaled chemical compositions and α-enhanced compositions with
[α/Fe] ≈ 0.4. The relation between the metal mass fraction ``Z`` and the logarithmic
metal abundance [M/H] is the same for the solar-scaled and α-enhanced models,
although the relation between the iron abundance and the overall metal abundance is
not the same -- [Fe/H] = [M/H] - 0.75 * [α/Fe]; see, e.g., Equation 4 of [Vazdekis2015](@citet).

The solar protostellar chemical mixture for these models was calibrated to
reproduce solar photospheric observations via a forward modeling approach
(see section 4 of [Pietrinferni2004](@citet)).
The distribution of heavy metals is taken from [Grevesse1993](@citet), which the authors
state is minimally different from [Grevesse1998](@citet).

```jldoctest
julia> using StellarTracks.BaSTIv1: BaSTIv1Chemistry, X, Y, Z, X_phot, Y_phot, Z_phot, MH;

julia> chem = BaSTIv1Chemistry();

julia> X(chem) + Y(chem) + Z(chem) ≈ 1 # solar protostellar values
true

julia> X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1 # solar photospheric values
true

julia> MH(chem, Z(chem) * 0.1) ≈ -0.9659287781921233
true

julia> Z(chem, -0.9659287781921233) ≈ Z(chem) * 0.1
true
```
"""
struct BaSTIv1Chemistry <: AbstractChemicalMixture end
X(mix::BaSTIv1Chemistry) = 1 - Y(mix) - Z(mix) # Protostellar abundance of calibrated solar model
X_phot(::BaSTIv1Chemistry) = 0.737993 # Calculated from Y_phot = 0.244, (Z/X)_phot = 0.0244
Y(::BaSTIv1Chemistry) = 0.2734 # Protostellar abundance of calibrated solar model
Y_phot(::BaSTIv1Chemistry) = 0.244 # Assumed present-day photospheric abundance
Y_p(::BaSTIv1Chemistry) = 0.245  # Second paragraph, section 5.1, Pietrinferni2004
Z(::BaSTIv1Chemistry) = 0.0198 # Protostellar abundance of calibrated solar model
Z_phot(::BaSTIv1Chemistry) = 0.018007 # Calculated from Y_phot = 0.244, (Z/X)_phot = 0.0244

Y(mix::BaSTIv1Chemistry, Zval) = Y_p(mix) +  14//10 * Zval # Second paragraph, section 5.1, Pietrinferni2004
# X generic
MH(mix::BaSTIv1Chemistry, Zval) = log10(Zval / X(mix, Zval)) - log10(Z_phot(mix) / X_phot(mix))
function Z(mix::BaSTIv1Chemistry, MHval)
    # Derivation in parsec code
    zoverx = exp10(MHval + log10(Z_phot(mix) / X_phot(mix)))
    γ = 14//10
    return (1 - Y_p(mix)) * zoverx / (1 + (1 + γ) * zoverx)
end

##########################################################################

# Data download, organization, etc.
include("init.jl")

##########################################################################

"""
    BaSTIv1Track(zval::Number, mass::Number, α_fe::Number, canonical::Bool, agb::Bool, η::Number)
`BaSTIv1Track` implements the [`AbstractTrack`](@ref StellarTracks.AbstractTrack)
interface for the older BaSTI stellar evolution library
[Pietrinferni2004,Pietrinferni2006,Pietrinferni2013](@citep).

Note that due to the organization of the BaSTIv1 data files, this method requires
constructing a [`BaSTIv1TrackSet`](@ref) and is therefore
not efficient if your aim is to construct multiple tracks with the same properties
but different masses. In this case, you should construct a [`BaSTIv1TrackSet`](@ref)
and call it with the masses you want, e.g.,
`ts = BaSTIv1TrackSet(0.0001, 0.0, true, true, 0.4); ts.([0.61, 0.82])`.

# Arguments
 - `zval::Number`: Metal mass fraction of stellar model.
 - `mass::Number`: Initial stellar mass of stellar model in solar masses.
 - `α_fe::Number`: [α/Fe] of stellar model.
 - `canonical::Bool`: Whether to use models with convective overshooting (`true`) or without (`false`).
 - `agb::Bool`: Whether to use models with the synthetic AGB extension (`true`) or without (`false`).
 - `η::Number`: Reimers mass loss parameter used to calculate the stellar model.

```jldoctest
julia> track = StellarTracks.BaSTIv1.BaSTIv1Track(1e-4, 1.2, 0.0, true, true, 0.4)
Canonical BaSTIv1Track with AGB extension, M_ini=1.2, MH=-2.2651979978981727, [α/Fe]=0.0, η=0.4, Z=0.0001, Y=0.24514, X=0.75476.
julia> track(9.0) # interpolate track at log10(age [yr]) = 9
(log_L = 0.703673172754993, log_Teff = 3.9489346126913416, log_g = 4.563808864709015)
```
"""
struct BaSTIv1Track{A,B,C} <: AbstractTrack
    data::Table{A}
    itp::B
    properties::C
end
# Constructor taking a subtable from one of the .jld2 reprocessed files
function BaSTIv1Track(data::Table, props)
    # Construct interpolator as a function of proper age
    itp = CubicSpline([SVector(values(i)[2:end]) for i in data],
                      deduplicate_knots!(data.star_age; move_knots=true);
                      cache_parameters=true)
    return BaSTIv1Track(data, itp, props)
end
# Constructor taking Z value, α_fe, canonical, initial stellar mass, loads Table, calls above method
function BaSTIv1Track(zval::Number, mass::Number, α_fe::Number, canonical::Bool, agb::Bool, η::Number)
    # For BaSTIv1, individual tracks are not saved, so we need to load a trackset
    props = (M = mass, Z = zval, α_fe = α_fe, canonical = canonical, agb = agb, η = η)
    data = BaSTIv1TrackSet(zval, α_fe, canonical, agb, η)(mass).data
    return BaSTIv1Track(data, props) # Method above
end
# Make Track callable with logAge to get log_L, log_Teff, and log_g as a NamedTuple
Base.keys(track::BaSTIv1Track) = columnnames(track.data)[2:end]
function (track::BaSTIv1Track)(logAge::Number)
    result = track.itp(exp10(logAge))
    return NamedTuple{keys(track)}(result)
end
Base.extrema(t::BaSTIv1Track) = log10.(extrema(t.itp.t))
mass(t::BaSTIv1Track) = t.properties.M
chemistry(::BaSTIv1Track) = BaSTIv1Chemistry()
MH(t::BaSTIv1Track) = MH(chemistry(t), Z(t))
Z(t::BaSTIv1Track) = t.properties.Z
Y(t::BaSTIv1Track) = Y(chemistry(t), Z(t))
X(t::BaSTIv1Track) = 1 - Y(t) - Z(t)
post_rgb(t::BaSTIv1Track) = length(t.data) > eep_idxs.HE_BEG
Base.eltype(t::BaSTIv1Track) = typeof(t.properties.Z)
function Base.show(io::IO, mime::MIME"text/plain", t::BaSTIv1Track)
    print(io, """$(ifelse(t.properties.canonical, "Canonical", "Non-canonical")) BaSTIv1Track $(ifelse(t.properties.agb, "with AGB extension", "without AGB extension")), M_ini=$(mass(t)), MH=$(MH(t)), [α/Fe]=$(t.properties.α_fe), η=$(t.properties.η), Z=$(Z(t)), Y=$(Y(t)), X=$(X(t)).""")
end

##########################################################################

"""
    BaSTIv1TrackSet(zval::Number, α_fe::Number=0, canonical::Bool=true, agb::Bool=false, η::Number=0.4)
`BaSTIv1TrackSet` implements the [`AbstractTrackSet`](@ref StellarTracks.AbstractTrackSet)
interface for the older BaSTI stellar evolution library
[Pietrinferni2004,Pietrinferni2006,Pietrinferni2013](@citep).

# Arguments
 - `zval::Number`: Metal mass fraction of stellar model.
 - `α_fe::Number = 0`: [α/Fe] of stellar model.
 - `canonical::Bool = true`: Whether to use models with convective overshooting (`true`) or without (`false`).
 - `agb::Bool = false`: Whether to use models with the synthetic AGB extension (`true`) or without (`false`).
 - `η::Number = 0.4`: Reimers mass loss parameter used to calculate the stellar model.

```jldoctest
julia> ts = StellarTracks.BaSTIv1.BaSTIv1TrackSet(1e-3, 0.0, true, false, 0.4)
Canonical BaSTIv1TrackSet without AGB extension, MH=-1.263953315381199, [α/Fe]=0.0, η=0.4, Z=0.001, Y=0.24640000006649643, 1999 EEPs and 40 initial stellar mass points.

julia> ts(1.01) # Interpolate track at new initial mass
Canonical BaSTIv1Track without AGB extension, M_ini=1.01, MH=-1.263953315381199, [α/Fe]=0.0, η=0.4, Z=0.001, Y=0.24640000006649643, X=0.7525999998860061.

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
# Given a set of properties, load the correct models and then call below method
function BaSTIv1TrackSet(zval::Number, α_fe::Number=0, canonical::Bool=true, agb::Bool=false, η::Number=0.4)
    # Validate arguments
    _validate_params(zval, α_fe, canonical, agb, η)
    # Convert to strings
    zval = _parse_zval(zval)
    α_fe = _parse_α_fe(α_fe)
    η = _parse_η(η)
    bfile = @datadep_str(joinpath("BaSTIv1", "basti_v1.jld2"))
    group = ifelse(canonical, "canonical", "noncanonical") * "/" *
            ifelse(agb, "AGB", "noAGB") * "/" *
            α_fe * "/" * zval * "/" * η
    data = JLD2.load(bfile, group)
    # data will now have whatever data types were originally saved into the jld2 file
    # We will promote to track_type here
    # data = Table(eep = data.eep, m_ini = convert(Vector{track_type}, data.m_ini),
    #              logAge = convert(Vector{track_type}, data.logAge),
    #              logL = convert(Vector{track_type}, data.logL),
    #              Teff = convert(Vector{track_type}, data.Teff),
    #              logg = convert(Vector{track_type}, data.logg))
    return BaSTIv1TrackSet(data, parse(track_type, zval), parse(track_type, α_fe), canonical, agb,
                           parse(track_type, η))
end
function BaSTIv1TrackSet(data::Table, zval::Number, α_fe::Number, canonical::Bool, agb::Bool, η::Number)
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

    Threads.@threads for i in eachindex(eeps)
    # for i in eachindex(eeps)
        eep = eeps[i]
        tmpdata = data[findall(Base.Fix1(==, eep), data.eep)] # Performance optimization
        # Sort by age, which will be the independent variable in the AMR interpolation
        idxs = sortperm(tmpdata.logAge)
        tmpdata = tmpdata[idxs]
        # Now keep only entries that define a unique mapping between age and m_ini
        tmpdata = tmpdata[uniqueidx(tmpdata.logAge)]

        # We'll start the interpolation at the most massive model and enforce a monotonic
        # decrease in stellar initial mass with increasing age at fixed EEP while EEP < eep_idxs[3]
        if eep < eep_idxs.MS_TO
            _, p1 = findmax(tmpdata.m_ini)
            idxs = p1:lastindex(tmpdata)
            goodidxs = diff(tmpdata.m_ini[idxs]) .< 0
            goodidxs = vcat(true, goodidxs) # add true for first element as well
            tmpdata = tmpdata[idxs[goodidxs]]
        end

        # PCHIPInterpolation is a type of CubicHermiteSpline
        amrs[i] = PCHIPInterpolation(tmpdata.m_ini, tmpdata.logAge)
        # Sort by initial stellar mass for defining the other interpolations
        idxs = sortperm(tmpdata.m_ini)
        tmpdata = tmpdata[idxs]
        # Now interpolate mass against logte, logl, logg
        logte[i] = PCHIPInterpolation(tmpdata.logTe, tmpdata.m_ini)
        logl[i] = PCHIPInterpolation(tmpdata.logL, tmpdata.m_ini)
        tmp_logg = log10.(surface_gravity.(tmpdata.mass, radius.(exp10.(tmpdata.logTe), tmpdata.logL)))
        logg[i] = PCHIPInterpolation(tmp_logg, tmpdata.m_ini)
    end
    return BaSTIv1TrackSet(eeps, amrs,
                           (log_L = logl, log_Teff = logte, log_g = logg),
                           (Z = zval, α_fe = α_fe, canonical = canonical,
                            agb = agb, η = η, masses = unique(data.m_ini)))
end
function (ts::BaSTIv1TrackSet)(M::Number)
    props = (M = M, Z = Z(ts), α_fe = ts.properties.α_fe, canonical = ts.properties.canonical,
             agb = ts.properties.agb, η = ts.properties.η)
    nt = _generic_trackset_interp(ts, M)
    table = Table(NamedTuple{(:star_age, keys(nt)[2:end]...)}(tuple(exp10.(nt.logAge), values(nt)[2:end]...)))
    return BaSTIv1Track(table, props)
end
mass(ts::BaSTIv1TrackSet) = ts.properties.masses
chemistry(::BaSTIv1TrackSet) = BaSTIv1Chemistry()
Z(ts::BaSTIv1TrackSet) = ts.properties.Z
MH(ts::BaSTIv1TrackSet) = MH(chemistry(ts), Z(ts))
Y(ts::BaSTIv1TrackSet) = Y(chemistry(ts), Z(ts))
X(ts::BaSTIv1TrackSet) = 1 - Y(ts) - Z(ts)
post_rgb(t::BaSTIv1TrackSet) = true
Base.eltype(ts::BaSTIv1TrackSet) = typeof(ts.properties.Z)
function Base.show(io::IO, mime::MIME"text/plain", ts::BaSTIv1TrackSet)
    print(io, """$(ifelse(ts.properties.canonical, "Canonical", "Non-canonical")) BaSTIv1TrackSet $(ifelse(ts.properties.agb, "with AGB extension", "without AGB extension")), MH=$(MH(ts)), [α/Fe]=$(ts.properties.α_fe), η=$(ts.properties.η), Z=$(Z(ts)), Y=$(Y(ts)), $(length(ts.AMRs)) EEPs and $(length(mass(ts))) initial stellar mass points.""")
end

function isochrone(ts::BaSTIv1TrackSet, logAge::Number)
    eeps = Vector{Int}(undef, 0)
    track_extrema = extrema(mass(ts))
    interp_masses = Vector{eltype(ts)}(undef, 0)
    logte = similar(interp_masses)
    logl = similar(interp_masses)
    logg = similar(interp_masses)
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
                # if length(logl) > 0 && i < eep_idxs[4] && logli < last(logl)
                # # if length(logl) > 0 && (i >= eep_idxs[3] && i <= eep_idxs[4]) && logli < last(logl)
                #     # li = lastindex(logl)
                #     # deleteat!(eeps, li)
                #     # deleteat!(interp_masses, li)
                #     # deleteat!(logte, li)
                #     # deleteat!(logl, li)
                #     # deleteat!(logg, li)
                #     continue
                # end
                push!(eeps, i)
                push!(interp_masses, imass)
                push!(logl, logli)
                push!(logte, ts.interps.log_Teff[i](imass))
                push!(logg, ts.interps.log_g[i](imass))
            end
        end
    end
    return (eep = eeps, m_ini = interp_masses, logTe = logte, Mbol = Mbol.(logl),
            logg = logg, logL = logl)
end

##########################################################################

"""
    BaSTIv1Library(α_fe::Number=0, canonical::Bool=true, agb::Bool=false, η::Number=0.4)
`BaSTIv1Library` implements the
[`AbstractTrackLibrary`](@ref StellarTracks.AbstractTrackLibrary)
interface for the older BaSTI stellar evolution models presented in
[Pietrinferni2004,Pietrinferni2006,Pietrinferni2013](@citet).

# Optional Arguments
 - `α_fe::Number = 0`: [α/Fe] of stellar model.
 - `canonical::Bool = false`: Whether to use models with convective overshooting (`true`) or without (`false`).
 - `agb::Bool = false`: Whether to use models with the synthetic AGB extension (`true`) or without (`false`).
 - `η::Number = 0.4`: Reimers mass loss parameter used to calculate the stellar model.

If you construct an instance as `p = BaSTIv1Library(0.0, true, false, 0.4)`, it is callable as
`p(mh::Number, M::Number)` which returns an [`InterpolatedTrack`](@ref StellarTracks.InterpolatedTrack)
that interpolates between tracks to a specific metallicity ([M/H]) and initial stellar mass (`M`).

This type also supports isochrone construction
(see [isochrone](@ref StellarTracks.isochrone(::StellarTracks.BaSTIv1.BaSTIv1Library, ::Number, ::Number))).

# Examples
```jldoctest
julia> p = BaSTIv1Library(0.0, true, false, 0.4)
Structure of interpolants for the older BaSTI library of canonical stellar tracks without AGB extension, [α/Fe]=0.0, η=0.4. Valid range of metallicities is (-3.265322279085031, 0.5088071641379128).

julia> isochrone(p, 10.05, -2.01) isa NamedTuple
true

julia> p(-2.05, 1.05)
InterpolatedTrack with M_ini=1.05, MH=-2.05, Z=0.0001641003350386593, Y=0.24522974046905413, X=0.7546061591959072.
```
"""
struct BaSTIv1Library{A,B} <: AbstractTrackLibrary
    ts::A   # Vector of `TrackSet`s
    properties::B
end
chemistry(::BaSTIv1Library) = BaSTIv1Chemistry()
Z(p::BaSTIv1Library) = p.properties.Z # Z.(chemistry(p), p.MH)
MH(p::BaSTIv1Library) = MH.(chemistry(p), Z(p))
Y(p::BaSTIv1Library) = Y.(chemistry(p), Z(p))
X(p::BaSTIv1Library) = 1 .- Y(p) .- Z(p)
post_rgb(::BaSTIv1Library) = true
Base.eltype(p::BaSTIv1Library) = typeof(first(Z(p)))
Base.Broadcast.broadcastable(p::BaSTIv1Library) = Ref(p)
function Base.show(io::IO, mime::MIME"text/plain", p::BaSTIv1Library)
    print(io, """Structure of interpolants for the older BaSTI library of $(ifelse(p.properties.canonical, "canonical", "non-canonical")) stellar tracks $(ifelse(p.properties.agb, "with AGB extension", "without AGB extension")), [α/Fe]=$(p.properties.α_fe), η=$(p.properties.η). Valid range of metallicities is $(extrema(MH(p))).""")
end
function BaSTIv1Library(α_fe::Number=0, canonical::Bool=true, agb::Bool=false, η::Number=0.4)
    # Only η=0.4 models with no AGB extension were run for first and last value of zgrid
    if η ≈ 0.4 && ~agb
        zg = zgrid
    else
        zg = zgrid[begin+1:end-1]
    end
    # Make vector of tracksets
    ts = [BaSTIv1TrackSet(z, α_fe, canonical, agb, η) for z in zg]
    return BaSTIv1Library(ts, (Z = zg, α_fe = α_fe, canonical = canonical, agb = agb, η = η))
end

# Below is a stub for documentation,
# calls down to generic method in StellarTracks.jl main file
"""
    isochrone(p::BaSTIv1Library, logAge::Number, mh::Number)
Interpolates properties of the stellar tracks in the library at the requested logarithmic age (`logAge = log10(age [yr])`) and logarithmic metallicity [M/H] = `mh`. Returns a `NamedTuple` containing the properties listed below:
 - `eep`: Equivalent evolutionary points
 - `m_ini`: Initial stellar masses, in units of solar masses.
 - `logTe`: Base-10 logarithm of the effective temperature [K] of the stellar model.
 - `Mbol`: Bolometric luminosity of the stellar model.
 - `logg`: Surface gravity of the stellar model.
"""
isochrone(p::BaSTIv1Library, logAge::Number, mh::Number)

export BaSTIv1Track, BaSTIv1TrackSet, BaSTIv1Library, BaSTIv1Chemistry # Unique module exports
export mass, chemistry, X, Y, Z, X_phot, Y_phot, Z_phot, MH, post_rgb, isochrone # Export generic API methods

end # module


# function plot_iso(iso)
#     fig,ax1 = plt.subplots()
#     ax1.plot(iso.logTe, iso.logL, marker="o", markersize=2, markeredgecolor="k")
#     # ax1.scatter([iso.logTe[iso.eep .== eep][1] for eep in values(eep_idxs)[2:end]],
#     # 	    [iso.logL[iso.eep .== eep][1] for eep in values(eep_idxs)[2:end]],
#     # 	    c="k")
#     ax1.scatter([iso.logTe[searchsortedfirst(iso.eep, eep)] for eep in values(eep_idxs)],
# 	        [iso.logL[searchsortedfirst(iso.eep, eep)] for eep in values(eep_idxs)],
# 	        c="k")
#     ax1.set_xlim([3.82, 3.65])
#     ax1.set_ylim([-0.5, 1.7])
# end

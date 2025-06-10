"""StellarTracks.BASTIv2 provides access to the updated BaSTI stellar tracks first presented in [Hidalgo2018](@citet)."""
module BaSTIv2

# imports from parent module
using ..StellarTracks: AbstractChemicalMixture, AbstractTrack, AbstractTrackSet,
                       AbstractTrackLibrary, uniqueidx, Mbol, _generic_trackset_interp,
                       radius, surface_gravity
import ..StellarTracks: X, Y, Z, X_phot, Y_phot, Z_phot, MH, chemistry, mass, post_rgb, isochrone
using ..StellarTracks.BaSTIv1: _parse_α_fe

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

# """Number of secondary EEP points per primary EEP point. [Hidalgo2018](@citet) Table 4, called key points (KPs). Key point 1 starts at age = 1000 yr."""
# const eep_lengths = (PMS_BEG = 20,  # age = 1000 yr to end of deuterium burning
#                      PMS_MINL = 60, # first min(L), or nuclear energy dominates total energy
#                      MS_BEG = 100,  # beginning of MS
#                      # intermediate age MS; min(Teff) for high-mass (HM); X_c = 0.3 for low-mass (LM)
#                      IAMS = 300,       
#                      # MS_KINK = 61,  # MS kink point; maximum Teff along MS
#                      MS_TO = 360,     # terminal MS; maximum Teff along MS (TO point)
#                      SGB = 420,       # Maximum log(L) for HM, X_c = 0 for LM
#                      RG_START = 490,  # min(L) for HM, RGB start for LM
#                      RG_BUMP = 10,     # beginning core He fusion
#                      HE_2 = 150,      # core He fusion Y_c = 0.55
#                      HE_3 = 100,      # core He fusion Y_c = 0.50
#                      HE_4 = 100,      # core He fusion Y_c = 0.40
#                      HE_5 = 80,       # core He fusion Y_c = 0.20
#                      HE_6 = 80,       # core He fusion Y_c = 0.10
#                      END_CHEB = 140,  # end of core He burning; core He fusion Y_c = 0.0
#                      TPAGB_BEG = 150) # L_{CNO} > L_{3α} during AGB phase

"""Indices into track tables of primary EEP points for BaSTIv2 given by [Hidalgo2018](@citet) Table 4, called key points (KPs)."""
const eep_idxs = (PMS_BEG = 1,    # age = 1000 yr
                  PMS_MINL = 20,  # end of deuterium burning
                  PMS_MID = 60,   # first min(L), or nuclear energy dominates total energy
                  MS_BEG = 100,   # beginning of MS
                  # intermediate age MS; min(Teff) for high-mass (HM); X_c = 0.3 for low-mass (LM)
                  IAMS = 300,
                  MS_TO = 360,    # terminal MS; maximum Teff along MS (TO point)
                  SGB = 420,      # max(L) for HM, X_c = 0 for LM
                  RG_START = 490, # min(L) for HM, RGB start for LM
                  RG_BUMP_MINL = 860, # min(L) along RGB bump
                  RG_BUMP_MAXL = 890, # max(L) along RGB bump
                  RG_TIP = 1290, # max(L) after core H-burning, RGB tip
                  HE_BEG = 1300, # beginning of core He-burning
                  HE_2 = 1450,      # core He fusion Y_c = 0.55
                  HE_3 = 1550,      # core He fusion Y_c = 0.50
                  HE_4 = 1650,      # core He fusion Y_c = 0.40
                  HE_5 = 1730,       # core He fusion Y_c = 0.20
                  HE_6 = 1810,       # core He fusion Y_c = 0.30
                  END_CHEB = 1950,   # end core helium burning, Y_c = 0.0
                  CNO_BEG = 2100)    # CNO cycle energy > helium burning energy

const eep_lengths = NamedTuple{keys(eep_idxs)[begin:end-1]}(eep_idxs[i+1] - eep_idxs[i] for i in eachindex(values(eep_idxs))[begin:end-1])

# Indices where the different phases begin
# """1-based indices giving the EEP point at which each EEP phase begins."""
# const eep_idxs = NamedTuple{keys(eep_lengths)}((1, (cumsum(values(eep_lengths)[begin:end-1]) .+ 1)...))
"""Data type to parse the BaSTIv2 tracks as."""
const track_type = Float32
"""Available [Fe/H] values in the BaSTIv2 stellar track grid."""
const feh_grid = track_type[-3.2, -2.5, -2.2, -1.9, -1.7, -1.55, -1.4, -1.3, -1.2, -1.05, -0.9, -0.7, -0.6, -0.4, -0.3, -0.2, -0.1, 0.06, 0.15, 0.3, 0.45] # -0.1 is listed as -0.08 in paper but -0.1 in the online grids...
"""Initial stellar masses for the stellar tracks in the BaSTIv2 grid. These are uniform for all metallicities and also for alpha-enhanced grids."""
const massgrid = track_type[0.1, 0.12, 0.15, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0]
"""α-element enrichment parameters [α/Fe] available for the BaSTIv2 grid."""
const αFegrid = track_type[-0.2, 0.0, 0.4]

function _validate_params(feh::Number, α_fe::Number, canonical::Bool, diffusion::Bool, yp::Number, η::Number)
    if ~any(Base.Fix1(isapprox, feh), BaSTIv2.feh_grid)
        throw(ArgumentError("Input [Fe/H] $feh invalid; available [Fe/H] values are $(feh_grid)."))
    end
    if α_fe ≈ 0
        if canonical
            if diffusion
                throw(ArgumentError("Models with diffusion are not available for canonical models with scaled-solar abundances. Set `diffusion=false`."))                
            else
                if ~(yp ≈ 0.247)
                    throw(ArgumentError("Canonical models with scaled-solar abundance patterns are only available for primordial helium abundance `yp = 0.247`."))
                end
                if ~(η ≈ 0)
                    throw(ArgumentError("Canonical models with scaled-solar abundance patterns are only available with Reimers mass loss parameter `η = 0`."))
                end
            end
        else
            if diffusion
                if ~(yp ≈ 0.247)
                    throw(ArgumentError("Non-canonical models with diffusion and scaled-solar abundance patterns are only available for primordial helium abundance `yp = 0.247`."))
                end
                if ~(η ≈ 0.3)
                    throw(ArgumentError("Non-canonical models with diffusion and scaled-solar abundance patterns are only available with Reimers mass loss parameter `η = 0.3`."))
                end
                if feh ≈ 0.45
                    throw(ArgumentError("Non-canonical models with diffusion and scaled-solar abundance patterns are only available for [Fe/H] values <= 0.3, you requested $feh."))
                end
            else
                if ~(yp ≈ 0.247)
                    throw(ArgumentError("Non-canonical models without diffusion with scaled-solar abundance patterns are only available for primordial helium abundance `yp = 0.247`."))
                end
                if ~any(η .≈ (0.0, 0.3))
                    throw(ArgumentError("Canonical models without diffusion with scaled-solar abundance patterns are only available with Reimers mass loss parameters `η = (0.0, 0.3)`."))
                end
            end
        end
    elseif α_fe ≈ 0.4
        if canonical
            throw(ArgumentError("Canonical models are not available for α-enhanced models with [α/Fe]=$α_fe. Set `canonical=false`."))
        else
            if diffusion
                
                if ~(η ≈ 0.3)
                    throw(ArgumentError("α-enhanced models with [α/Fe]=$α_fe are only available for Reimers mass loss parameter `η=0.3`."))
                end
                if yp ≈ 0.247
                    if ~(feh < 0.06 || feh ≈ 0.06)
                        throw(ArgumentError("α-enhanced models with [α/Fe]=$α_fe and Y_p=$yp are only available for iron abundances up to [Fe/H] = 0.06, you requested $feh."))
                    end
                elseif yp ≈ 0.275
                    if ~(feh < -0.6 || feh ≈ -0.6)
                        throw(ArgumentError("α-enhanced models with [α/Fe]=$α_fe and Y_p=$yp are only available for iron abundances up to [Fe/H] = -0.6, you requested $feh."))
                    end                    
                elseif yp ≈ 0.3
                    if ~(feh < -0.1 || feh ≈ -0.1)
                        throw(ArgumentError("α-enhanced models with [α/Fe]=$α_fe and Y_p=$yp are only available for iron abundances up to [Fe/H] = -0.08, you requested $feh."))
                    end
                elseif yp ≈ 0.32
                    if ~(feh ≈ 0.06)
                        throw(ArgumentError("α-enhanced models with [α/Fe]=$α_fe and Y_p=$yp are only available for [Fe/H] = 0.06, you requested $feh."))
                    end
                end
                
            else
                throw(ArgumentError("Models without diffusion are not available for α-enhanced models with [α/Fe]=$α_fe. Set `diffusion=true`."))
            end
        end

    elseif α_fe ≈ -0.2

    else
        throw(ArgumentError("Provided [α/Fe] abundance $α_fe invalid; available values are (-0.2, 0.0, 0.4)."))
    end
    return true
end

"""
    _parse_feh(fehval::Number)
Convert a logarithmic abundance into a properly formatted string for use in loading the JLD2 file.
"""
_parse_feh(fehval::Number) = @sprintf("%1.2f", fehval)
"""
    _parse_yp(yp::Number)
Convert the primordial helium abundance into a properly formatted string for use in loading the JLD2 file.
"""
_parse_yp(yp::Number) = @sprintf("%1.3f", yp)
"""
    _parse_η(η::Number)
Convert the Reimers mass loss parameter into a properly formatted string for use in loading the JLD2 file.
"""
_parse_η(η::Number) = @sprintf("%1.1f", η)

##########################################################################
# BaSTIv2 chemistry
"""
    BaSTIv2Chemistry(α_fe::Number, yp::Number)
Returns a struct representing the BaSTIv2 chemical mixture model with [α/Fe] = `α_fe`
and primordial helium abundance `yp`. 
These BaSTI models, presented in [Hidalgo2018,Pietrinferni2021,Salaris2022,Pietrinferni2024](@citet),
include solar-scaled chemical compositions, α-enhanced compositions with
[α/Fe] ≈ 0.4, and α-depleted compositions with [α/Fe] ≈ -0.2.

The solar protostellar chemical mixture for these models was calibrated to
reproduce solar photospheric observations via a forward modeling approach
(see section 3 of [Hidalgo2018](@citet)). Most solar photospheric abundances
are taken from [Caffau2011](@citet).

As the BaSTIv2 grid was run with uniform [Fe/H] values but differing [α/Fe], the
metal mass fraction ``Z`` and the logarithmic metal abundance [M/H] are not uniform
for every [α/Fe] in the grid. We therefore need to know [α/Fe] so that we can convert
the [Fe/H] values to [M/H] and ``Z``. The α-enhanced library also includes models
with enhanced primordial helium abundance (`yp` here), so we require that information
as well.

For the [α/Fe] = 0.4 models presented in [Pietrinferni2021(@citet) and the [α/Fe] = -0.2 models
presented in [Pietrinferni2024](@citet), the α elements
O, Ne, Mg, Si, S, Ca, and Ti have all been uniformly modified with respect to the Fe abundance
relative to the [Caffau2011](@citet) heavy element distribution. As C and N are not changed,
[M/H] = [Fe/H] + 0.75 * [α/Fe]; see, e.g., Equation 4 of [Vazdekis2015](@citet).

```jldoctest
julia> using StellarTracks.BaSTIv2: BaSTIv2Chemistry, X, Y, Z, X_phot, Y_phot, Z_phot, MH;

julia> chem = BaSTIv2Chemistry(0.0, 0.247);

julia> X(chem) + Y(chem) + Z(chem) ≈ 1 # solar protostellar values
true

julia> X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1 # solar photospheric values
true

julia> isapprox(MH(chem, Z_phot(chem) * 0.1), -1; atol=0.01)
true

julia> isapprox(Z(chem, -1), Z_phot(chem) * 0.1; rtol=0.02)
true
```
"""
struct BaSTIv2Chemistry{T} <: AbstractChemicalMixture
    α_fe::T
    yp::T
end
BaSTIv2Chemistry(α_fe, yp) = BaSTIv2Chemistry(promote(α_fe, yp)...)
X(::BaSTIv2Chemistry) = 0.7133 # Protostellar abundance of calibrated solar model
X_phot(::BaSTIv2Chemistry) = 0.7362 # Calculated from Z_phot, Y_phot
Y(::BaSTIv2Chemistry) = 0.2695 # Protostellar abundance of calibrated solar model
Y_phot(::BaSTIv2Chemistry) = 0.2485 # Present-day photospheric abundance; Basu 2004
Y_p(mix::BaSTIv2Chemistry) = mix.yp # 0.247  # First paragraph, section 4, Hidalgo2018
Z(::BaSTIv2Chemistry) = 0.0172 # Protostellar abundance of calibrated solar model
Z_phot(::BaSTIv2Chemistry) = 0.0153 # Caffau2011

Y(mix::BaSTIv2Chemistry, Zval) = Y_p(mix) + 1.31 * Zval # First paragraph, section 4, Hidalgo2018
# X generic

# Its just [M/H] = [Fe/H] + 0.75 * [α/Fe]
MH(mix::BaSTIv2Chemistry, Zval) = log10(Zval / X(mix, Zval)) - log10(Z_phot(mix) / X_phot(mix))
function Z(mix::BaSTIv2Chemistry, MHval)
    # Derivation in parsec code
    zoverx = exp10(MHval + log10(Z_phot(mix) / X_phot(mix)))
    γ = 1.31
    return (1 - Y_p(mix)) * zoverx / (1 + (1 + γ) * zoverx)
end

##########################################################################

# Data download, organization, etc.
include("init.jl")

##########################################################################
"""
    BaSTIv2Track(feh::Number, mass::Number, α_fe::Number, canonical::Bool, diffusion::Bool, yp::Number, η::Number)
`BaSTIv2Track` implements the [`AbstractTrack`](@ref StellarTracks.AbstractTrack)
interface for the updated BaSTI stellar evolution library
[Hidalgo2018,Pietrinferni2021,Salaris2022,Pietrinferni2024](@citep).

Note that due to the organization of the BaSTIv2 data files, this method requires
constructing a [`BaSTIv2TrackSet`](@ref) and is therefore
not efficient if your aim is to construct multiple tracks with the same properties
but different masses. In this case, you should construct a [`BaSTIv2TrackSet`](@ref)
and call it with the masses you want, e.g.,
`ts = BaSTIv2TrackSet(-2.2, 0.0, true, true, 0.247, 0.3); ts.([0.61, 0.82])`.

# Arguments
 - `feh::Number`: [Fe/H] of stellar model
 - `mass::Number`: Initial stellar mass of stellar model in solar masses
 - `α_fe::Number`: [α/Fe] of stellar model
 - `canonical::Bool`: Whether to use models with convective overshooting (`true`) or without (`false`).
 - `diffusion::Bool`: Whether to use models with atomic diffusion (`true`) or without (`false`).
 - `yp::Number`: Primordial helium abundance assumed for stellar model.
 - `η::Number`: Reimers mass loss parameter used to calculate the stellar model.

Note that this function takes the input metallicity as [Fe/H], which is not equal to [M/H]
when considering models with [α/Fe] != 0 `α_fe != 0`.
```jldoctest
julia> track = StellarTracks.BaSTIv2.BaSTIv2Track(-2.2, 0.81, 0.0, false, true, 0.247, 0.3)
Non-canonical BaSTIv2Track with diffusion, M_ini=0.81, [M/H]=-2.2, [Fe/H]=-2.2, [α/Fe]=0.0, Z=9.870952533235687e-5, Y=0.24712930947818537, X=0.7527719809964823, Y_p=0.247, η=0.3.

julia> track(9.0) # interpolate track at log10(age [yr]) = 9
(log_L = -0.14188867667533753, log_Teff = 3.799111830003377, log_g = 4.639490166221381)
```
"""
struct BaSTIv2Track{A,B,C} <: AbstractTrack
    data::Table{A}
    itp::B
    properties::C
end
# Constructor taking a subtable from one of the .jld2 reprocessed files
function BaSTIv2Track(data::Table, props)
    # Construct interpolator as a function of proper age
    itp = CubicSpline([SVector(values(i)[2:end]) for i in data],
                      deduplicate_knots!(data.star_age; move_knots=true))
    return BaSTIv2Track(data, itp, props)
end
# Constructor taking Z value, α_fe, canonical, diffusion, initial stellar mass, loads Table, calls above method
function BaSTIv2Track(feh::Number, mass::Number, α_fe::Number, canonical::Bool, diffusion::Bool,
                      yp::Number, η::Number)
    _validate_params(feh, α_fe, canonical, diffusion, yp, η)
    # For BaSTIv2, individual tracks are not saved, so we need to load a trackset
    props = (M = mass, feh = feh, α_fe = α_fe, canonical = canonical, diffusion = diffusion, yp = yp, η = η)
    data = BaSTIv2TrackSet(feh, α_fe, canonical, diffusion, yp, η)(mass).data
    return BaSTIv2Track(data, props) # Method above
end
# Make Track callable with logAge to get log_L, log_Teff as a NamedTuple
function (track::BaSTIv2Track)(logAge::Number)
    result = track.itp(exp10(logAge))
    return NamedTuple{columnnames(track.data)[2:end]}(result)
end
Base.extrema(t::BaSTIv2Track) = log10.(extrema(t.itp.t))
mass(t::BaSTIv2Track) = t.properties.M
chemistry(t::BaSTIv2Track) = BaSTIv2Chemistry(t.properties.α_fe, t.properties.yp)
MH(t::BaSTIv2Track) = t.properties.feh + 3//4 * t.properties.α_fe # MH(chemistry(t), Z(t))
Z(t::BaSTIv2Track) = Z(chemistry(t), MH(t)) # t.properties.Z
Y(t::BaSTIv2Track) = Y(chemistry(t), Z(t))
X(t::BaSTIv2Track) = 1 - Y(t) - Z(t)
post_rgb(t::BaSTIv2Track) = length(t.data) > eep_idxs.HE_BEG
Base.eltype(t::BaSTIv2Track) = typeof(t.properties.feh)
function Base.show(io::IO, mime::MIME"text/plain", t::BaSTIv2Track)
    print(io, """$(ifelse(t.properties.canonical, "Canonical", "Non-canonical")) BaSTIv2Track $(ifelse(t.properties.diffusion, "with diffusion", "without diffusion")), M_ini=$(mass(t)), [M/H]=$(MH(t)), [Fe/H]=$(t.properties.feh), [α/Fe]=$(t.properties.α_fe), Z=$(Z(t)), Y=$(Y(t)), X=$(X(t)), Y_p=$(Y_p(chemistry(t))), η=$(t.properties.η).""")
end

##########################################################################

"""
    BaSTIv2TrackSet(feh::Number, α_fe::Number=0, canonical::Bool=false,
                    diffusion::Bool=true, yp::Number=0.247, η::Number=0.3)
`BaSTIv2TrackSet` implements the [`AbstractTrackSet`](@ref StellarTracks.AbstractTrackSet)
interface for the updated BaSTI stellar evolution library
[Hidalgo2018,Pietrinferni2021,Salaris2022,Pietrinferni2024](@citep).

# Arguments
 - `feh::Number`: [Fe/H] of stellar model

# Optional Arguments
 - `α_fe::Number = 0`: [α/Fe] of stellar model
 - `canonical::Bool = false`: Whether to use models with convective overshooting (`true`) or without (`false`).
 - `diffusion::Bool = true`: Whether to use models with atomic diffusion (`true`) or without (`false`).
 - `yp::Number = 0.247`: Primordial helium abundance assumed for stellar model.
 - `η::Number = 0.3`: Reimers mass loss parameter used to calculate the stellar model.

```jldoctest
julia> ts = StellarTracks.BaSTIv2.BaSTIv2TrackSet(-2.2, 0.0, false, true, 0.247, 0.3)
Non-canonical BaSTIv2TrackSet with diffusion, [M/H]=-1.9, [Fe/H]=-1.9, [α/Fe]=0.0, Z=0.00019689206760795184, Y=0.2472579286085664, Y_p=0.247, η=0.3, 2099 EEPs and 56 initial stellar mass points.

julia> ts(1.01) # Interpolate track at new initial mass
Non-canonical BaSTIv2Track with diffusion, M_ini=1.01, [M/H]=-1.9, [Fe/H]=-1.9, [α/Fe]=0.0, Z=0.00019689206760795184, Y=0.2472579286085664, X=0.7525451793238256, Y_p=0.247, η=0.3.

julia> isochrone(ts, 10.0) isa NamedTuple # Interpolate isochrone at `log10(age [yr]) = 10`
true
```
"""
struct BaSTIv2TrackSet{A <: AbstractVector{<:Integer},
                       B <: AbstractVector{<:AbstractInterpolation},
                       C,
                       D} <: AbstractTrackSet
    eeps::A # EEP points; vector of integers
    AMRs::B # age-mass relations; vector of interpolators, same length as eeps
    interps::C
    properties::D
end
# Given a metallicity, α-enhancement, and canonical status, load the correct models and then call below method
function BaSTIv2TrackSet(feh::Number, α_fe::Number=0, canonical::Bool=false, diffusion::Bool=true,
                         yp::Number=0.247, η::Number=0.3)
    # Validate arguments
    _validate_params(feh, α_fe, canonical, diffusion, yp, η)
    # Convert to strings
    feh = _parse_feh(feh)
    α_fe = _parse_α_fe(α_fe)
    yp = _parse_yp(yp)
    η = _parse_η(η)
    dd_path = @datadep_str("BaSTIv2")
    bfile = joinpath(dd_path, "basti_v2.jld2")
    group = ifelse(canonical, "canonical", "noncanonical") * "/" *
            ifelse(diffusion, "diffusion", "nodiffusion") * "/" *
            α_fe * "/" * feh * "/" * yp * "/" * η
    println(group)
    data = JLD2.load(bfile, group)
    # data will now have whatever data types were originally saved into the jld2 file
    # We will promote to track_type here
    # data = Table(eep = data.eep, m_ini = convert(Vector{track_type}, data.m_ini),
    #              logAge = convert(Vector{track_type}, data.logAge),
    #              logL = convert(Vector{track_type}, data.logL),
    #              Teff = convert(Vector{track_type}, data.Teff),
    #              logg = convert(Vector{track_type}, data.logg))
    return BaSTIv2TrackSet(data, parse(track_type, feh), parse(track_type, α_fe), canonical, diffusion,
                           parse(track_type, yp), parse(track_type, η))
end
function BaSTIv2TrackSet(data::Table, feh::Number, α_fe::Number, canonical::Bool, diffusion::Bool,
                         yp::Number, η::Number)
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
        _, p1 = findmax(tmpdata.m_ini)
        idxs = p1:lastindex(tmpdata)
        if eep < eep_idxs.MS_TO
            goodidxs = diff(tmpdata.m_ini[idxs]) .< 0
            goodidxs = vcat(true, goodidxs) # add true for first element as well
        else
            goodidxs = trues(length(idxs))
        end
        tmpdata = tmpdata[idxs[goodidxs]]

        # PCHIPInterpolation is a type of CubicHermiteSpline
        amrs[i] = PCHIPInterpolation(tmpdata.m_ini, tmpdata.logAge)
        # Sort by initial stellar mass for defining the other interpolations
        idxs = sortperm(tmpdata.m_ini)
        tmpdata = tmpdata[idxs]
        # Now interpolate mass against logte, logl, logg
        logte[i] = PCHIPInterpolation(tmpdata.logTe, tmpdata.m_ini)
        logl[i] = PCHIPInterpolation(tmpdata.logL, tmpdata.m_ini)
        # Calculate logg from mass, temperature, luminosity
        logg[i] = PCHIPInterpolation(log10.(surface_gravity.(tmpdata.mass, radius.(exp10.(tmpdata.logTe), tmpdata.logL))),
                                     tmpdata.m_ini)
    end
    return BaSTIv2TrackSet(eeps, amrs,
                           (log_L = logl, log_Teff = logte, log_g = logg),
                           (feh = feh, α_fe = α_fe, canonical = canonical, diffusion = diffusion,
                            yp = yp, η = η, masses = unique(data.m_ini)))
end
function (ts::BaSTIv2TrackSet)(M::Number)
    props = (M = M, feh = MH(ts), α_fe = ts.properties.α_fe, canonical = ts.properties.canonical,
             diffusion=ts.properties.diffusion, yp=ts.properties.yp, η=ts.properties.η)
    nt = _generic_trackset_interp(ts, M)
    table = Table(NamedTuple{(:star_age, keys(nt)[2:end]...)}(tuple(exp10.(nt.logAge), values(nt)[2:end]...)))
    return BaSTIv2Track(table, props)
end
mass(ts::BaSTIv2TrackSet) = ts.properties.masses
chemistry(ts::BaSTIv2TrackSet) = BaSTIv2Chemistry(ts.properties.α_fe, ts.properties.yp)
Z(ts::BaSTIv2TrackSet) = Z(chemistry(ts), MH(ts)) # ts.properties.Z
MH(ts::BaSTIv2TrackSet) = ts.properties.feh + 3//4 * ts.properties.α_fe # MH(chemistry(ts), Z(ts))
Y(ts::BaSTIv2TrackSet) = Y(chemistry(ts), Z(ts))
X(ts::BaSTIv2TrackSet) = 1 - Y(ts) - Z(ts)
post_rgb(t::BaSTIv2TrackSet) = true
Base.eltype(ts::BaSTIv2TrackSet) = typeof(ts.properties.feh)
function Base.show(io::IO, mime::MIME"text/plain", ts::BaSTIv2TrackSet)
    print(io, """$(ifelse(ts.properties.canonical, "Canonical", "Non-canonical")) BaSTIv2TrackSet $(ifelse(ts.properties.diffusion, "with diffusion", "without diffusion")), [M/H]=$(MH(ts)), [Fe/H]=$(ts.properties.feh), [α/Fe]=$(ts.properties.α_fe), Z=$(Z(ts)), Y=$(Y(ts)), Y_p=$(Y_p(chemistry(ts))), η=$(ts.properties.η), $(length(ts.AMRs)) EEPs and $(length(mass(ts))) initial stellar mass points.""")
end

function isochrone(ts::BaSTIv2TrackSet, logAge::Number)
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
    BaSTIv2Library(α_fe::Number=0, canonical::Bool=false, diffusion::Bool=true,
                   yp::Number=0.247, η::Number=0.3)
`BaSTIv2Library` implements the
[`AbstractTrackLibrary`](@ref StellarTracks.AbstractTrackLibrary)
interface for the updated BaSTI stellar evolution models presented in
[Hidalgo2018,Pietrinferni2021,Salaris2022,Pietrinferni2024](@citep).

# Optional Arguments
 - `α_fe::Number = 0`: [α/Fe] of stellar model
 - `canonical::Bool = false`: Whether to use models with convective overshooting (`true`) or without (`false`).
 - `diffusion::Bool = true`: Whether to use models with atomic diffusion (`true`) or without (`false`).
 - `yp::Number = 0.247`: Primordial helium abundance assumed for stellar model.
 - `η::Number = 0.3`: Reimers mass loss parameter used to calculate the stellar model.

If you construct an instance as `p = BaSTIv2Library(0.0, false)`, it is callable as
 - `p(mh::Number)` to interpolate the full library to a new metallicity
   (returning a [`BaSTIv2TrackSet`](@ref)), or
 - `p(mh::Number, M::Number)` to interpolate the tracks to a specific metallicity
   and initial stellar mass (returning a [`BaSTIv2Track`](@ref)).

This type also supports isochrone construction
(see [isochrone](@ref StellarTracks.isochrone(::StellarTracks.BaSTIv2.BaSTIv2Library, ::Number, ::Number))).

```jldoctest
julia> p = BaSTIv2Library(0.0, false, true, 0.247, 0.3)
Structure of interpolants for the updated BaSTI library of non-canonical stellar tracks with diffusion, [α/Fe]=0.0, Y_p=0.247, η=0.3. Valid range of metallicities is [Fe/H] = $(extrema(feh_grid)), [M/H] = $(extrema(feh_grid)).

julia> isochrone(p, 10.05, -2.01) isa NamedTuple
true
```
"""
struct BaSTIv2Library{A,B} <: AbstractTrackLibrary
    ts::A   # Vector of `TrackSet`s
    properties::B
    # feh::B    # Vector of [Fe/H] for each TrackSet
    # α_fe::C # α-enhancement
    # canonical::Bool
    # diffusion::Bool
    # yp::D
end
# Interpolation to get a TrackSet with metallicity MH
function (ts::BaSTIv2Library)(mh::Number)
    error("Not yet implemented.")
end
# Interpolation to get a Track with mass M and metallicity MH
function (ts::BaSTIv2Library)(mh::Number, M::Number)
    error("Not yet implemented.")
end
chemistry(p::BaSTIv2Library) = BaSTIv2Chemistry(p.properties.α_fe, p.properties.yp)
Z(p::BaSTIv2Library) = Z.(chemistry(p), MH(p))
MH(p::BaSTIv2Library) = p.properties.feh .+ 3//4 .* p.properties.α_fe # MH.(chemistry(p), Z(p))
Y(p::BaSTIv2Library) = Y.(chemistry(p), Z(p))
X(p::BaSTIv2Library) = 1 .- Y(p) .- Z(p)
post_rgb(::BaSTIv2Library) = true
Base.eltype(p::BaSTIv2Library) = typeof(first(MH(p)))
Base.Broadcast.broadcastable(p::BaSTIv2Library) = Ref(p)
function Base.show(io::IO, mime::MIME"text/plain", p::BaSTIv2Library)
    print(io, """Structure of interpolants for the updated BaSTI library of $(ifelse(p.properties.canonical, "canonical", "non-canonical")) stellar tracks $(ifelse(p.properties.diffusion, "with diffusion", "without diffusion")), [α/Fe]=$(p.properties.α_fe), Y_p=$(p.properties.yp), η=$(p.properties.η). Valid range of metallicities is [Fe/H] = $(extrema(p.properties.feh)), [M/H] = $(extrema(MH(p))).""")
end
function BaSTIv2Library(α_fe::Number=0, canonical::Bool=false, diffusion::Bool=true,
                        yp::Number=0.247, η::Number=0.3)
    # Not all values in feh_grid are supported for every combination of parameters
    # Need try/catch to gracefully handle unsupported feh values
    ts = [begin
              try
                  BaSTIv2TrackSet(feh, α_fe, canonical, diffusion, yp, η)
              catch
              end
          end for feh in feh_grid]
    if any(isnothing, ts)
        ts = filter(!isnothing, ts)
    end
    return BaSTIv2Library(ts, (feh = [tt.properties.feh for tt in ts],
                               α_fe = convert(track_type, α_fe),
                               canonical = canonical, diffusion = diffusion,
                               yp = convert(track_type, yp), η = convert(track_type, η)))
end

# Below is a stub for documentation,
# calls down to generic method in StellarTracks.jl main file
"""
    isochrone(p::BaSTIv2Library, logAge::Number, mh::Number)
Interpolates properties of the stellar tracks in the library at the requested logarithmic age (`logAge = log10(age [yr])`) and logarithmic metallicity [M/H] = `mh`. Returns a `NamedTuple` containing the properties listed below:
 - `eep`: Equivalent evolutionary points
 - `m_ini`: Initial stellar masses, in units of solar masses.
 - `logTe`: Base-10 logarithm of the effective temperature [K] of the stellar model.
 - `Mbol`: Bolometric luminosity of the stellar model.
 - `logg`: Surface gravity of the stellar model.
"""
isochrone(p::BaSTIv2Library, logAge::Number, mh::Number)

export BaSTIv2Track, BaSTIv2TrackSet, BaSTIv2Library, BaSTIv2Chemistry # Unique module exports
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

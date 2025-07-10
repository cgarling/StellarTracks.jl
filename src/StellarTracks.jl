module StellarTracks

using ArgCheck: @argcheck
using TypedTables: Table

# For BCs.jl
using BolometricCorrections: AbstractBCTable, MISTBCGrid, PHOENIXYBCGrid, filternames
import BolometricCorrections: AbstractChemicalMixture, X, X_phot, Y, Y_phot, Z, Z_phot, Y_p, MH, chemistry
using DataInterpolations: PCHIPInterpolation
import Tables

# Top-level API definitions
# Track
""" Abstract supertype representing individual stellar tracks.
Different concrete implementations are used for different
stellar libraries, but all subtypes of `AbstractTrack` are guaranteed to support
the generic API described in the documentation.

The main way to get the properties of the
modeled star when it has a particular age, assuming `track` is a valid concrete instance,
is to call the track with the logarithmic age `track(log10(age))` where `age` is in years.
This will return a `NamedTuple` containing the available properties (e.g., logL, logTe)
of the modeled star at that age. Not all libraries will have the same properties available. """
abstract type AbstractTrack end 
Base.Broadcast.broadcastable(t::AbstractTrack) = Ref(t)
(track::AbstractTrack)(logAge::AbstractArray{<:Number}) = Table(track(la) for la in logAge)
# function (track::AbstractTrack)(logAge::AbstractArray{<:Number})
#     # result = track.itp(logAge)
#     # return Table(NamedTuple{(:logTe, :Mbol, :logg)}(i) for i in result)
#     # return Table(NamedTuple{(:logTe, :Mbol, :logg)}(track.itp(la)) for la in logAge)
#     # Conversion to Table is slightly slow ~40ns 
#     return Table(track(la) for la in logAge)
# end
"""
    Base.extrema(t::AbstractTrack)
Returns the minimum and maximum logarithmic age (`log10(age [yr])`) of the stellar model.
"""
function Base.extrema(t::AbstractTrack) end
"""
    Base.keys(t::AbstractTrack)
Returns a tuple of symbols that will be used to set the keys of the `NamedTuple`s returned by `t(logAge::Number)`.
"""
function Base.keys(t::AbstractTrack) end
"""
    mass(t::AbstractTrack)
Returns the initial stellar mass of the modeled star in solar masses. """
function mass(t::AbstractTrack) end
"""
    chemistry(t::AbstractTrack)
Returns an instance of a subtype of
[`AbstractChemicalMixture`](@extref BolometricCorrections.AbstractChemicalMixture)
valid for the provided track `t` that can be used in other methods.
"""
function chemistry(t::AbstractTrack) end
"""
    X(t::AbstractTrack)
Returns the hydrogen mass fraction of the modeled star. """
function X(t::AbstractTrack) end
"""
    Y(t::AbstractTrack)
Returns the helium mass fraction of the modeled star. """
function Y(t::AbstractTrack) end
"""
    Z(t::AbstractTrack)
Returns the metal mass fraction of the modeled star. """
function Z(t::AbstractTrack) end
"""
    MH(t::AbstractTrack)
Returns the logarithmic metal abundance of the modeled star, defined as [M/H] = log(Z/X) - log(Z⊙/X⊙). """
function MH(t::AbstractTrack) end
"""
    post_rgb(t::AbstractTrack)
Returns `true` if the track includes post-RGB evolution, `false` otherwise. """
function post_rgb(t::AbstractTrack) end
###############################################
# TrackSet
""" Abstract supertype representing a set of stellar tracks
with common properties and supports
common operations like interpolating between tracks and
calculating isochrones. Specifically, it is assumed all individual
tracks in a track set have identical initial metallicity and chemical
composition. Subtypes of `AbstractTrackSet`
are guaranteed to support the generic API described in the documentation.

Concrete instances are callable with an initial stellar mass (in solar masses),
returning an interpolated track at the requested mass. """
abstract type AbstractTrackSet end
Base.Broadcast.broadcastable(ts::AbstractTrackSet) = Ref(ts)
# Generic function to interpolate a trackset to a new initial stellar mass
# This returns a NamedTuple; concrete subtypes should define methods to
# construct the correct AbstractTrack type from this NamedTuple.
function _generic_trackset_interp(ts::AbstractTrackSet, M::Number)
    # Validate that mass is in range
    # throw(DomainError(M, "Requested mass $M is outside the valid range $(extrema(mass(ts))) for the track set."))
    m_min, m_max = extrema(mass(ts))
    @argcheck m_min <= M <= m_max "Requested mass $M is outside the valid range $(extrema(mass(ts))) for the track set."
    interps = ts.interps
    results = Vector{Vector{eltype(ts)}}(undef, length(interps)) # +1 for age column
    # Find EEP points where the requested mass is valid
    good_idx = findall(ii -> begin
                                  ee = extrema(first(interps)[ii])
                                  (M >= ee[1]) && (M <= ee[2])
                             end, eachindex(first(interps)))
    # Loop over unique values that are being interpolated (logg, Teff, etc), one interp per property
    for i in eachindex(values(interps))
        interps_i = interps[i]
        results[i] = [interps_i[j](M) for j in good_idx]
    end
    # Nearly all computation time is spent here -- faster interpolation inversion
    # or *maybe* root-finding on existing interpolation would speed this up
    ages = [begin
                sortidx = sortperm(amr.u)
                PCHIPInterpolation(amr.t[sortidx], amr.u[sortidx])(M)
            end for amr in ts.AMRs[good_idx]]
    sortidx = sortperm(ages)
    return NamedTuple{(:logAge, keys(interps)...)}(tuple(ages[sortidx], (r[sortidx] for r in results)...))
end
"""
    mass(ts::AbstractTrackSet)
Returns the initial stellar masses (in solar masses) of the individual tracks contained in the track set. """
function mass(ts::AbstractTrackSet) end
"""
    chemistry(ts::AbstractTrack)
Returns an instance of a subtype of
[`AbstractChemicalMixture`](@extref BolometricCorrections.AbstractChemicalMixture)
valid for the provided track set `ts` that can be used in other methods.
"""
function chemistry(ts::AbstractTrackSet) end
"""
    X(ts::AbstractTrackSet)
Returns the common hydrogen mass fraction of the tracks contained in the track set. """
function X(ts::AbstractTrackSet) end
"""
    Y(ts::AbstractTrackSet)
Returns the common helium mass fraction of the tracks contained in the track set. """
function Y(ts::AbstractTrackSet) end
"""
    Z(ts::AbstractTrackSet)
Returns the common metal mass fraction of the tracks contained in the track set. """
function Z(ts::AbstractTrackSet) end
"""
    MH(ts::AbstractTrackSet)
Returns the common logarithmic metal abundance of the tracks contained in the track set, defined as [M/H] = log(Z/X) - log(Z⊙/X⊙). """
function MH(ts::AbstractTrackSet) end
"""
    post_rgb(ts::AbstractTrackSet)
Returns `true` if the tracks in the track set include post-RGB evolution, `false` otherwise. """
function post_rgb(ts::AbstractTrackSet) end
"""
    isochrone(ts::AbstractTrackSet, logAge::Number)
Interpolates properties of the stellar tracks in the track set at
the requested logarithmic age (`logAge = log10(age [yr])`).
Returns a `NamedTuple` containing the properties; different libraries
may have different properties available. If `result = isochrone(...)`,
EEP points are generally available as `result.eep` and the corresponding
initial stellar masses are `result.m_ini`.

```jldoctest
julia> using StellarTracks.MIST

julia> ts = MISTTrackSet(-1, 0); # Load set of MIST tracks with [M/H] = -1, vvcrit=0

julia> isochrone(ts, 10.0) isa NamedTuple
true
```
"""
function isochrone(ts::AbstractTrackSet, logAge::Number) end

###############################################
# TrackLibrary
""" Abstract supertype for loading the *full* set of all
stellar tracks available from a given library (e.g., PARSEC).
For some libraries (e.g., MIST), there are different model sets
for other stellar input parameters (e.g., rotation `vvcrit` in
MIST) that are not interpolated over -- these should be provided
as arguments to the concrete type constructors.

These are typically assembled as a collection of individual
track sets, with one track set constructed for each unique
set of stellar chemical compositions (i.e., [M/H]).
For more details, see the documentation for
[`AbstractTrackSet`](@ref StellarTracks.AbstractTrackSet).

Concrete subtypes must implement an isochrone interpolation method.
A generic method for the call signature
[`isochrone(tracklib::AbstractTrackLibrary, logAge::Number, mh::Number)`](@ref isochrone(::AbstractTrackLibrary, ::Number, ::Number))
is provided. This generic method should work as long as the concrete
subtype follows the [chemistry API](@ref chemistry_api) and has a field `tracklib.ts`
containing the individual tracksets that make up the library. The tracksets must support
isochrone interpolation via a method `isochrone(ts::<your_type_here>, logAge::Number)`.

Concrete instances must support interpolation of stellar tracks
at user-provided metallicities and initial stellar masses. A generic method
is implemented via the signature
`(tracklib::AbstractTrackLibrary)(mh::Number, M::Number)`. This method
interpolates the `tracklib` at logarithmic metallicity [M/H] = `mh` and stellar
initial mass `M`, returning an [`InterpolatedTrack`](@ref StellarTracks.InterpolatedTrack)
that can be called to evaluate the track."""
abstract type AbstractTrackLibrary end

"""
    chemistry(tl::AbstractTrackLibrary)
Returns an instance of a subtype of
[`AbstractChemicalMixture`](@extref BolometricCorrections.AbstractChemicalMixture)
valid for the provided track library `tl` that can be used in other methods.
"""
function chemistry(tl::AbstractTrackLibrary) end
"""
    X(tl::AbstractTrackLibrary)
Returns the hydrogen mass fractions of the track sets 
contained in the track library. """
function X(tl::AbstractTrackLibrary) end
"""
    Y(tl::AbstractTrackLibrary)
Returns the helium mass fractions of the track sets
contained in the track library. """
function Y(tl::AbstractTrackLibrary) end
"""
    Z(tl::AbstractTrackLibrary)
Returns the metal mass fractions of the track sets
contained in the track library. """
function Z(tl::AbstractTrackLibrary) end
"""
    MH(tl::AbstractTrackLibrary)
Returns the logarithmic metal abundances of the track sets 
contained in the track library, 
defined as [M/H] = log(Z/X) - log(Z⊙/X⊙). """
function MH(tl::AbstractTrackLibrary) end
"""
    post_rgb(tl::AbstractTrackLibrary)
Returns `true` if any of the track sets
in the track library include post-RGB evolution, 
`false` otherwise. """
function post_rgb(tl::AbstractTrackLibrary) end

"""
    InterpolatedTrack(track0, track1, track1_prefac, track2_prefac) <: AbstractTrack

Type allowing for linear interpolation between two stellar tracks of identical 
mass but different metallicity (`track0` and `track1`). Simple linear 
interpolation between the two tracks is applied. The track can be evaluated 
at a particular logarithmic age by calling `(track::InterpolatedTrack)(logAge::Number)`. 

The constructor for this type is considered internal and users should not construct this type directly.
Rather, tracks of this type are constructed from track libraries by calling them with 
the signature `(tracklib::AbstractTrackLibrary)(mh::Number, M::Number)` where 
`mh` is the logarithmic metallicity [M/H] and `M` is the initial stellar mass in solar masses.

```jldoctest
julia> using StellarTracks.PARSEC: PARSECLibrary

julia> p = PARSECLibrary();

julia> track = p(-2.05, 1.05)
InterpolatedTrack with M_ini=1.05, MH=-2.05, Z=0.00013856708164357998, Y=0.24874664940532557, X=0.7511147835130308.

julia> track(9.0)
(logTe = 3.8820487347062302, Mbol = 3.7411721770340987, logg = 4.521853108813156, C_O = 0.0)
```
"""
struct InterpolatedTrack{A <: AbstractTrack, B <: Number} <: AbstractTrack
    track0::A
    track1::A
    track0_prefac::B
    track1_prefac::B
end
function (track::InterpolatedTrack)(logAge::Number)
    # Just let it fail if logAge errors for one of the two tracks
    result0 = track.track0(logAge)
    result1 = track.track1(logAge)
    result = values(result0) .* track.track0_prefac .+ values(result1) .* track.track1_prefac
    return NamedTuple{keys(track.track0)}(result)
end
function Base.extrema(track::InterpolatedTrack)
    ext0 = extrema(track.track0)
    ext1 = extrema(track.track1)
    return (max(ext0[1], ext1[1]), min(ext0[2], ext1[2]))
end
Base.keys(track::InterpolatedTrack) = keys(track.track0)
Base.eltype(track::InterpolatedTrack) = promote_type(eltype(track.track0), eltype(track.track1))
mass(track::InterpolatedTrack) = mass(track.track0) == mass(track.track1) ? mass(track.track0) : mass(track.track0) * track.track0_prefac + mass(track.track1) * track.track1_prefac
chemistry(track::InterpolatedTrack) = chemistry(track.track0)
MH(track::InterpolatedTrack) = MH(track.track0) * track.track0_prefac + MH(track.track1) * track.track1_prefac
Z(track::InterpolatedTrack) = Z(chemistry(track), MH(track))
Y(track::InterpolatedTrack) = Y(chemistry(track), Z(track))
X(track::InterpolatedTrack) = X(chemistry(track), Z(track))
post_rgb(track::InterpolatedTrack) = all(post_rgb, (track.track0, track.track1)) ? true : false
function Base.show(io::IO, mime::MIME"text/plain", t::InterpolatedTrack{A}) where A
    print(io, "InterpolatedTrack with M_ini=$(mass(t)), MH=$(MH(t)), Z=$(Z(t)), Y=$(Y(t)), X=$(X(t)).")
end
# function Base.show(io::IO, mime::MIME"text/plain", t::InterpolatedTrack{A}) where A
#     print(io, "Interpolation between two `$(split(string(A), "{")[1])`s with M_ini=$(mass(t)), MH=$(MH(t)), Z=$(Z(t)), Y=$(Y(t)), X=$(X(t)).")
# end

# Make AbstractTrackLibrary callable to construct interpolated track
function (tracklib::AbstractTrackLibrary)(mh::Number, M::Number)
    ts = tracklib.ts # vector of tracksets that make up tracklib
    mhvec = MH(tracklib)
    min_mh, max_mh = extrema(mhvec)
    if mh < min_mh || mh > max_mh
        throw(DomainError(mh, "Requested metallicity [M/H]=$mh is outside the valid range for the stellar track library of $(extrema(mhvec))."))
    end

    if !issorted(mhvec)
        idxs = sortperm(mhvec)
        ts = view(ts, idxs)
        mhvec = view(mhvec, idxs)
    end
    
    idx = searchsortedfirst(mhvec, mh)
    if mhvec[idx] ≈ mh # Requested mh matched in grid
        return ts[idx](M) # Assume a method (trackset)(M) exists
    else
        # Need to interpolate between two tracks with same mass, different metallicity
        track0 = ts[idx-1](M)
        track1 = ts[idx](M)
        mh0 = mhvec[idx-1]
        mh1 = mhvec[idx]
        return InterpolatedTrack(track0, track1, ((mh1-mh)/(mh1-mh0)), ((mh-mh0)/(mh1-mh0)))
    end
end

# This generic method was originally developed for MIST;
# should work for PARSEC as well after refactoring
"""
    isochrone(tracklib::AbstractTrackLibrary, logAge::Number, mh::Number)
Interpolates properties of the stellar tracks in the track library
at the requested logarithmic age (`logAge = log10(age [yr])`) and
metallicity [M/H] = `mh`. Returns a `NamedTuple` containing the
properties; different libraries may have different properties available.
If `result = isochrone(...)`, EEP points are generally available as
`result.eep` and the corresponding initial stellar masses are `result.m_ini`.

```jldoctest
julia> using StellarTracks.MIST

julia> p = MISTLibrary(0.0); # Load the library of MIST tracks with vvcrit=0

julia> isochrone(p, 10.0, -1.65) isa NamedTuple
true
```
"""
function isochrone(tracklib::AbstractTrackLibrary, logAge::Number, mh::Number)
    ts = tracklib.ts
    mhvec = MH(tracklib)
    min_mh, max_mh = extrema(mhvec)
    if mh < min_mh || mh > max_mh
        throw(DomainError(mh, "Requested metallicity [M/H]=$mh is outside the valid range for the stellar track library of $(extrema(mhvec))."))
    end

    if !issorted(mhvec)
        idxs = sortperm(mhvec)
        ts = view(ts, idxs)
        mhvec = view(mhvec, idxs)
    end

    # searchsortedfirst returns the index of the first value in mhvec greater than or
    # equivalent to mh. If mh is greater than all values in mhvec, returns lastindex(mhvec) + 1.
    # We have already checked bounds so we know min_mh < mh < max_mh
    idx = searchsortedfirst(mhvec, mh)
    if mhvec[idx] ≈ mh # Requested mh matched in grid
        return isochrone(ts[idx], logAge)
    end

    # mh is valid, so need to interpolate isochrone as a function of [M/H]
    # According to Marigo2017, the interpolations (at least for BCs) in Z or [M/H]
    # are linear, so the BCs overall are computed by a
    # 3D linear interpolation in logTe x logg x [M/H] space.
    # VandenBerg2014 suggests linear interpolation is sufficient for Z or [M/H] interpolation,
    # although cubic interpolation was used in VandenBerg2012.
    # Most pre-computed grid seems more uniform in [M/H] than they are in Z, so I think it might
    # be a good idea to do linear interpolation in [M/H].
    
    # Evaluate isochrones on either side of intermediate point
    y0 = isochrone(ts[idx-1], logAge)
    y1 = isochrone(ts[idx], logAge)
    if length(first(y0)) == 0 || length(first(y1)) == 0
        throw(DomainError(logAge, "No valid EEPs were found for the requested logarithmic age $logAge."))
    end
    # Get indices into y0 and y1 that correspond to the overlapping EEP points
    y0_idxs = Vector{Int}(undef, 0)
    y1_idxs = similar(y0_idxs)
    good_eeps = similar(y0_idxs)
    for (y0_idx, eep) in enumerate(y0.eep)
        y1_idx = searchsortedfirst(y1.eep, eep)
        if y1_idx < lastindex(y1.eep) + 1 # Match found
            push!(y0_idxs, y0_idx)
            push!(y1_idxs, y1_idx)
            push!(good_eeps, eep)
        end
    end
    # Get isochrone keys, removing EEP since that is fixed
    goodkeys = filter(Base.Fix1(!==, :eep), keys(y0))
    # Perform linear interpolation in _interp_kernel to establish a function
    # barrier to improve performance since some of the types of the variables
    # aren't known at runtime
    result = NamedTuple{goodkeys}(_interp_kernel(goodkeys, y0, y1, idx, y0_idxs, y1_idxs, mh, mhvec))
    # Concatenate interpolated result with valid EEP points
    return (eep = good_eeps, result...)
end
# This does linear interpolation for isochrone; function barrier for type stability
_interp_kernel(goodkeys, y0, y1, idx, y0_idxs, y1_idxs, x, xvec) =
    ((y0[key][y0_idxs] .* (xvec[idx] - x) .+ y1[key][y1_idxs] .* (x - xvec[idx-1])) ./ (xvec[idx] - xvec[idx-1]) for key in goodkeys)

###############################################
# Utilities
"""
    uniqueidx(v) = unique(Base.Fix1(getindex, v), eachindex(v))
Returns the indices of the first occurrences of unique elements in the array `v`. """
uniqueidx(v) = unique(Base.Fix1(getindex, v), eachindex(v)) # utility
"""
    Mbol(logL::Number, solmbol::Number=4.74)
Returns the bolometric magnitude corresponding to the provided logarithmic
bolometric luminosity `logL` which is provided in units of solar luminosities
(e.g., `logL = log10(L / L⊙)`). This is given by `Mbol⊙ - 2.5 * logL`; the zeropoint
of bolometric magnitude scale is defined by the solar bolometric magnitude, which you
can specify as the second argument. The default (4.74) was recommended by
IAU [Resolution B2](@cite Mamajek2015).
"""
Mbol(logL::Number, solmbol::Number=474//100) = solmbol - 5 * logL / 2
"""
    logL(Mbol::Number, solmbol::Number=4.74)
Returns the logarithmic bolometric luminosity in units of solar luminosities
(e.g., `logL = log10(L / L⊙)`) corresponding to the provided bolometric
magnitude. This is given by `(Mbol⊙ - Mbol) / 2.5`; the zeropoint
of bolometric magnitude scale is defined by the solar bolometric magnitude, 
which you can specify as the second argument. The default (4.74) was recommended
by IAU [Resolution B2](@cite Mamajek2015).
"""
logL(Mbol::Number, solmbol::Number=474//100) = (solmbol - Mbol) / 5 * 2

"""
    radius(Teff::Number, logl::Number)
Returns the radius of a star in units of solar radii
given its effective temperature `Teff` in Kelvin
and the logarithm of its luminosity in units of solar luminosities
(e.g., `logl = log10(L / L⊙)`).

Assumes solar properties following [IAU 2015 Resolution B3](@cite Mamajek2015a).
```jldoctest
julia> isapprox(StellarTracks.radius(5772, 0.0), 1.0; rtol=0.001) # For solar Teff and logL, radius ≈ 1
true
```
"""
radius(Teff::Number, logl::Number) = sqrt(exp10(logl) / 4 / π / (Teff^2)^2) * 11810222860206199 // 100000000
# radius(Teff::Number, logl::Number) = sqrt(exp10(logl) / 4 / π / Teff^4) * 1.1810222860206199e8
# radius(Teff, logl) = sqrt(exp10(logl) * 3.828e26 / 4 / π / 5.6703744191844294e-8 / Teff^4) / 6.957e8
# radius(Teff, logl) = sqrt(exp10(logl) * UnitfulAstro.Lsun / 4 / π / PhysicalConstants.CODATA2022.StefanBoltzmannConstant / (Teff * UnitfulAstro.K)^4) |> UnitfulAstro.m
"""
    surface_gravity(M, R)
Returns the surface gravity of a star in cgs units `cm / s^2`
given its mass in solar masses and radius in solar radii.

Assumes solar properties following [IAU 2015 Resolution B3](@cite Mamajek2015a).
```jldoctest
julia> isapprox(StellarTracks.surface_gravity(1, 1), 27420; rtol=0.001) # For solar M and R, g ≈ 27430 cm / s^2
true
```
"""
surface_gravity(M::Number, R::Number) = M / R^2 * 27420011165737313 // 1000000000000
# surface_gravity(M, R) = 27420.011165737313 * M / R^2
# surface_gravity(M, R) = PhysicalConstants.CODATA2022.G * M * UnitfulAstro.Msun / (R * UnitfulAstro.Rsun)^2 |> UnitfulAstro.cm / UnitfulAstro.s^2

###############################################
# Include files containing submodules for different track libraries
include(joinpath("parsec", "parsec.jl"))
using .PARSEC

include(joinpath("mist", "mist.jl"))
using .MIST

include(joinpath("basti", "v1", "basti_v1.jl"))
using .BaSTIv1

include(joinpath("basti", "v2", "basti_v2.jl"))
using .BaSTIv2

# Include bolometric correction-related functionality
include("BCs.jl")

# Common API exports
export mass, chemistry, X, Y, Z, MH, post_rgb, isochrone
# Submodule exports
export PARSECLibrary, MISTLibrary, BaSTIv1Library, BaSTIv2Library

end

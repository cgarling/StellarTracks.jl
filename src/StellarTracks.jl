module StellarTracks

using ArgCheck: @argcheck
using TypedTables: Table

# For BCs.jl
using BolometricCorrections: AbstractBCTable, MISTBCGrid, filternames
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
"""
    Base.extrema(t::AbstractTrack)
Returns the minimum and maximum logarithmic age (`log10(age [yr])`) of the stellar model.
"""
function Base.extrema(t::AbstractTrack) end
"""
    mass(t::AbstractTrack)
Returns the initial stellar mass of the modeled star in solar masses. """
function mass(t::AbstractTrack) end
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
    interp_keys = keys(interps)
    interp_length = length(first(interps))
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
Subtypes of `AbstractTrackLibrary` are guaranteed to support
the generic API described in the documentation.

Concrete instances can be called to interpolate tracks
at other metallicities and initial stellar masses --
specifics of this behavior depend on the library. """
abstract type AbstractTrackLibrary end

# This generic method was originally developed for MIST;
# should work for PARSEC as well after refactoring
"""
    isochrone(tl::AbstractTrackLibrary, logAge::Number, mh::Number)
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
function isochrone(p::AbstractTrackLibrary, logAge::Number, mh::Number)
    mh_vec = MH(p)
    idx = findfirst(Base.Fix1(≈, mh), mh_vec) # Will be === nothing if no entry in MH(p) is ≈ mh
    # If input mh is represented in base grid, no mh interpolation needed
    if !isnothing(idx) # idx !== nothing
        return isochrone(p.ts[idx], logAge)
    end
    # Check mh is in valid range
    min_mh, max_mh = extrema(mh_vec)
    if mh < min_mh || mh > max_mh
        throw(DomainError(mh, "Requested metallicity [M/H]=$mh is outside the valid range for MIST library of $(extrema(mh_vec))."))
    end
    # mh is valid, so need to interpolate isochrone as a function of [M/H]
    # According to Marigo2017, the interpolations (at least for BCs) in Z or [M/H]
    # are linear, so the BCs overall are computed by a
    # 3D linear interpolation in logTe x logg x [M/H] space.
    # VandenBerg2014 suggests linear interpolation is sufficient for Z or [M/H] interpolation,
    # although cubic interpolation was used in VandenBerg2012.
    # Most pre-computed grid seems more uniform in [M/H] than they are in Z, so I think it might
    # be a good idea to do linear interpolation in [M/H].
    
    # searchsortedfirst returns the index of the first value in mh_vec greater than or
    # equivalent to mh. If mh is greater than all values in mh_vec, returns lastindex(mh_vec) + 1.
    # We have already checked bounds so we know min_mh < mh < max_mh
    idx = searchsortedfirst(mh_vec, mh)
    # Evaluate isochrones on either side of intermediate point
    y0 = isochrone(p.ts[idx-1], logAge)
    y1 = isochrone(p.ts[idx], logAge)
    if length(first(y0)) == 0 || length(first(y1)) == 0
        throw(DomainError(logAge, "No valid EEPs were found for the requested logarithmic age $logAge."))
    end
    # Get intersection of valid EEPs from each isochrone
    min_eep = max(first(y0.eep), first(y1.eep))
    max_eep = min(last(y0.eep), last(y1.eep))
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
    result = NamedTuple{goodkeys}(_interp_kernel(goodkeys, y0, y1, idx, y0_idxs, y1_idxs, mh, mh_vec))
    # Concatenate interpolated result with valid EEP points
    return (eep = good_eeps, result...)
end
# This does linear interpolation for isochrone
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
Mbol(logL::Number, solmbol::Number=4.74) = solmbol - 5 * logL / 2
"""
    logL(Mbol::Number, solmbol::Number=4.74)
Returns the logarithmic bolometric luminosity in units of solar luminosities
(e.g., `logL = log10(L / L⊙)`) corresponding to the provided bolometric
magnitude. This is given by `(Mbol⊙ - Mbol) / 2.5`; the zeropoint
of bolometric magnitude scale is defined by the solar bolometric magnitude, 
which you can specify as the second argument. The default (4.74) was recommended
by IAU [Resolution B2](@cite Mamajek2015).
"""
logL(Mbol::Number, solmbol::Number=4.74) = (solmbol - Mbol) / 5 * 2

###############################################
# Include files containing submodules for different track libraries
include(joinpath("parsec", "parsec.jl"))
using .PARSEC

include(joinpath("mist", "mist.jl"))
using .MIST

include(joinpath("basti", "v1", "basti_v1.jl"))
using .BaSTIv1

# Include bolometric correction-related functionality
include("BCs.jl")

# Common API exports
export mass, chemistry, X, Y, Z, MH, post_rgb, isochrone
# Submodule exports
export PARSECLibrary, MISTLibrary, BaSTIv1Library

end

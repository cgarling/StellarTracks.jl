module StellarTracks

using TypedTables: Table

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
Interpolates properties of the stellar tracks in the track set at the requested logarithmic age (`logAge = log10(age [yr])`). Returns a `NamedTuple` containing the properties; different libraries may have different properties available. If `result = isochrone(...)`, EEP points are generally available as `result.eep` and the corresponding initial stellar masses are `result.m_ini`.
"""
function isochrone(ts::AbstractTrackSet, logAge::Number) end

###############################################
# TrackLibrary
""" Abstract supertype for loading the *full* set of all
stellar tracks available from a given library (e.g., PARSEC).
These are typically assembled as a collection of individual
track sets, with one track set constructed for each unique
set of stellar chemical compositions.
For more details, see the documentation for
 [`AbstractTrackSet`](@ref StellarTracks.AbstractTrackSet). 
Subtypes of `AbstractTrackLibrary` are guaranteed to support
the generic API described in the documentation.

Concrete instances can be called to interpolate tracks
at other metallicities and initial stellar masses --
specifics of this behavior depend on the library. """
abstract type AbstractTrackLibrary end


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
include("parsec/parsec.jl")

# Common API exports
export mass, X, Y, Z, MH, post_rgb, isochrone
# Submodule exports
export PARSECLibrary

end

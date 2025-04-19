```@meta
CurrentModule = StellarTracks
```

# [Overview](@id overview)

## Individual Tracks
An individual stellar track, containing the time evolution of properties for a single star, is represented by our [`AbstractTrack`](@ref StellarTracks.AbstractTrack) type. Different concrete implementations are used for different stellar evolution libraries, but all support the common access API defined below.

### API
```@docs
StellarTracks.AbstractTrack
mass(::StellarTracks.AbstractTrack)
X(::StellarTracks.AbstractTrack)
Y(::StellarTracks.AbstractTrack)
Z(::StellarTracks.AbstractTrack)
MH(::StellarTracks.AbstractTrack)
post_rgb(::StellarTracks.AbstractTrack)
```

### Concrete Implementations
 - [`PARSEC.PARSECTrack`](@ref)

## Track Sets
We define a "track set" to be a set of individual stellar tracks that share common properties (typically initial metallicity). By grouping these tracks together, we can interpolate between tracks in the set, create isochrones, and perform other similar operations.

### API
```@docs
StellarTracks.AbstractTrackSet
mass(::StellarTracks.AbstractTrackSet)
X(::StellarTracks.AbstractTrackSet)
Y(::StellarTracks.AbstractTrackSet)
Z(::StellarTracks.AbstractTrackSet)
MH(::StellarTracks.AbstractTrackSet)
post_rgb(::StellarTracks.AbstractTrackSet)
isochrone(::StellarTracks.AbstractTrackSet, ::Number)
```

### Concrete Implementations
 - [`PARSEC.PARSECTrackSet`](@ref)

## Track Libraries
For ease of use, our main entry point is the [`AbstractTrackLibrary`](@ref StellarTracks.AbstractTrackLibrary), which loads and organizes all stellar tracks available from a given library (e.g., PARSEC). Once all tracks have been loaded, subsets with common chemical compositions can be extracted. Individual tracks and isochrones can also be interpolated directly from the library instance. However, the interfaces for these interpolations are not generic as not all libraries offer the same variations in initial chemistry -- e.g., most will offer variation in total metallicity (i.e., ``Z``), but some also include variation in ``\alpha``-element abundances. These methods are documented separately for each stellar library under their unique pages in the left panel.

### API
```@docs
StellarTracks.AbstractTrackLibrary
```

### Concrete Implementations
 - [`PARSEC.PARSECLibrary`](@ref)

## Utilities
```@docs
StellarTracks.Mbol
StellarTracks.logL
```

## Index
```@index
```

```@meta
CurrentModule = StellarTracks
```

# [API](@id api)

## [Chemistry](@id chemistry_api)
We utilize the [chemistry API from BolometricCorrections.jl](@extref BolometricCorrections chemistry_api) to describe the chemistry of the stellar track libraries.

## Individual Tracks
An individual stellar track, containing the time evolution of properties for a single star, is represented by our [`AbstractTrack`](@ref StellarTracks.AbstractTrack) type. Different concrete implementations are used for different stellar evolution libraries, but all support the common access API defined below.

### API
```@docs
StellarTracks.AbstractTrack
gridname(::StellarTracks.AbstractTrack)
Base.extrema(::StellarTracks.AbstractTrack)
Base.keys(::StellarTracks.AbstractTrack)
mass(::StellarTracks.AbstractTrack)
chemistry(::StellarTracks.AbstractTrack)
X(::StellarTracks.AbstractTrack)
Y(::StellarTracks.AbstractTrack)
Z(::StellarTracks.AbstractTrack)
MH(::StellarTracks.AbstractTrack)
FeH(::StellarTracks.AbstractTrack)
alphaFe(::StellarTracks.AbstractTrack)
post_rgb(::StellarTracks.AbstractTrack)
```

### Concrete Implementations
 - [`PARSEC.PARSECTrack`](@ref)
 - [`MIST.MISTv1Track`](@ref)
 - [`BaSTIv1.BaSTIv1Track`](@ref)
 - [`BaSTIv2.BaSTIv2Track`](@ref)
 - [`StellarTracks.InterpolatedTrack`](@ref)

```@docs
StellarTracks.InterpolatedTrack
```

## Track Sets
We define a "track set" to be a set of individual stellar tracks that share common properties (typically initial metallicity). By grouping these tracks together, we can interpolate between tracks in the set, create isochrones, and perform other similar operations.

For each EEP, age–mass relations are constructed by scanning tracks from highest to lowest initial mass and retaining only ages strictly greater than the last retained age. Equal ages retain the higher-mass sample. This top-down selection follows the canonical approach in [Dotter2016](@citet): physical nonmonotonic regions are represented by a single monotonic selection rather than multiple masses per EEP. Stellar properties use the same selected tracks. This approximation can change EEP mass intervals and therefore IMF-weighted phase counts. Selection occurs in memory and does not modify the saved track files. If selection removes a mass from every EEP of a track set, that mass is excluded from the supported mass grid rather than advertising a track with no valid points.

At each EEP, the logarithm of initial mass is interpolated linearly in logarithmic age between the selected tracks. Because the interpolation weights depend only on the track masses, the initial masses of an isochrone cannot decrease between neighboring EEPs that use the same selected tracks. Cubic interpolants do not have this property.

### API
```@docs
StellarTracks.AbstractTrackSet
gridname(::StellarTracks.AbstractTrackSet)
mass(::StellarTracks.AbstractTrackSet)
chemistry(::StellarTracks.AbstractTrackSet)
X(::StellarTracks.AbstractTrackSet)
Y(::StellarTracks.AbstractTrackSet)
Z(::StellarTracks.AbstractTrackSet)
MH(::StellarTracks.AbstractTrackSet)
FeH(::StellarTracks.AbstractTrackSet)
alphaFe(::StellarTracks.AbstractTrackSet)
post_rgb(::StellarTracks.AbstractTrackSet)
isochrone(::StellarTracks.AbstractTrackSet, ::Number)
isochrone(::StellarTracks.AbstractTrackSet, ::StellarTracks.AbstractBCTable, ::Number)
```

### Concrete Implementations
 - [`PARSEC.PARSECTrackSet`](@ref)
 - [`MIST.MISTv1TrackSet`](@ref)
 - [`BaSTIv1.BaSTIv1TrackSet`](@ref)
 - [`BaSTIv2.BaSTIv2TrackSet`](@ref)

## Track Libraries
For ease of use, our main entry point is the [`AbstractTrackLibrary`](@ref StellarTracks.AbstractTrackLibrary), which loads and organizes all stellar tracks available from a given library (e.g., PARSEC). Once all tracks have been loaded, subsets with common chemical compositions can be extracted. Individual tracks and isochrones can also be interpolated directly from the library instance. However, the interfaces for these interpolations are not generic as not all libraries offer the same variations in initial chemistry -- e.g., most will offer variation in total metallicity (i.e., ``Z``), but some also include variation in ``\alpha``-element abundances. These methods are documented separately for each stellar library under their unique pages in the left panel.

### API
```@docs
StellarTracks.AbstractTrackLibrary
gridname(::StellarTracks.AbstractTrackLibrary)
chemistry(::StellarTracks.AbstractTrackLibrary)
X(::StellarTracks.AbstractTrackLibrary)
Y(::StellarTracks.AbstractTrackLibrary)
Z(::StellarTracks.AbstractTrackLibrary)
MH(::StellarTracks.AbstractTrackLibrary)
FeH(::StellarTracks.AbstractTrackLibrary)
alphaFe(::StellarTracks.AbstractTrackLibrary)
post_rgb(::StellarTracks.AbstractTrackLibrary)
isochrone(::StellarTracks.AbstractTrackLibrary, ::Number, ::Number)
isochrone(::StellarTracks.AbstractTrackLibrary, ::StellarTracks.AbstractBCTable, ::Number, ::Number)
```

### Concrete Implementations
 - [`PARSEC.PARSECLibrary`](@ref)
 - [`MIST.MISTv1Library`](@ref)
 - [`BaSTIv1.BaSTIv1Library`](@ref)
 - [`BaSTIv2.BaSTIv2Library`](@ref)
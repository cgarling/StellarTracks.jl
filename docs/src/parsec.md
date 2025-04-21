```@meta
CurrentModule = StellarTracks
ShareDefaultModule = true
```

# [PARSEC](@id PARSEC)

Here we describe the interface we provide to the PARSEC v1.2S library of stellar evolutionary tracks. PARSEC specific code is housed in the `PARSEC` submodule, which can be accessed as

```julia
using StellarTracks.PARSEC # load all exported methods
using StellarTracks.PARSEC: PARSECLibrary, X, Y, Z # load specific methods
```

The main paper describing the PARSEC family of stellar models is [Bressan2012](@citet), but the library of stellar models has been expanded over the years to add and improve coverage of various parameter spaces. A non-exhaustive list of papers presenting the PARSEC models up to V1.2S is provided below.
 - [Bressan2012](@citet) is the first paper presenting the PARSEC models.
 - [Chen2014](@citet) presented an improved calibration for low-mass stars.
 - [Tang2014](@citet) and [Chen2015](@citet) presented new models of high-mass stars from 14 to 350 M⊙ with metallicities from ``0.001 \le Z \le 0.004``.
 - [Rosenfield2016](@citet) formulated equivalent evolutionary points (EEPs) for the PARSEC models to support use in isochrone interpolation routines. **Our implementation uses their data products.**
 - [Marigo2017](@citet) augmented the PARSEC models with COLIBRI models of the thermally pulsating asymptotic giant branch phase (TP-AGB).
 - [Pastorelli2019,Pastorelli2020](@citet) compared the COLIBRI TP-AGB models to observations of the SMC and LMC, respectively.

The full list of relevant papers maintained by the group is available [here](https://ui.adsabs.harvard.edu/public-libraries/jSpa1621SGW2mMpPRaRP4w).

As we use the PARSEC V1.2S tracks augmented with EEP points by [Rosenfield2016](@citet), we do not currently support the more recent PARSEC V2.0 tracks [Nguyen2022,Costa2025](@citep). We hope to add V2.0 in the future, but doing so would require new measurements of the EEP points which we do not presently support.

## Data Acquisition

This package handles downloading and pre-processing of the EEP tracks produced by [Rosenfield2016](@citet) (available [here](https://github.com/philrosenfield/padova_tracks)) using [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl). The main access point we provide is [`PARSECLibrary`](@ref StellarTracks.PARSEC.PARSECLibrary), which will load and make available the full library of stellar tracks. The first time you call this method, you will be prompted to download the required data files. The total data volume is ~66 MB after processing. Information on customizing the install location [here](https://www.oxinabox.net/DataDeps.jl/stable/z10-for-end-users/). The data can be uninstalled by running `using DataDeps; rm(datadep"PARSECv1.2S"; recursive=true)`. With all the tracks available, we are able to perform operations like interpolating isochrones at any age and metallicity with the PARSEC parameter space.

## Examples
Load the full PARSEC library, which is downloaded via DataDeps.jl if not already available.
```@example
using StellarTracks.PARSEC
p = PARSECLibrary()
```

Use the [`PARSEC.PARSECLibrary`](@ref) to interpolate an isochrone, at `log10(age [yr]) = 10.05` and metal mass fraction ``Z=0.001654``. The isochrone is returned as a `NamedTuple`.
```@example
isochrone(p, 10.05, 0.001654)
```

The `NamedTuple` returned by `isochrone()` can be converted to table types, like `TypedTables.Table` to simplify further use.
```@example
using TypedTables: Table
Table(isochrone(p, 10.05, 0.001654))
```

## Full Library
```@docs
StellarTracks.PARSEC.PARSECLibrary
isochrone(::StellarTracks.PARSEC.PARSECLibrary, ::Number, ::Number)
```

The full library is principally a set of [`PARSECTrackSet`](@ref StellarTracks.PARSEC.PARSECTrackSet)s, with one track set per unique chemical composition. All PARSEC models have scaled-solar chemical compositions, so they vary only in total metallicity (i.e., ``Z``). 

## Track Sets
```@docs
StellarTracks.PARSEC.PARSECTrackSet
```

Each track set is, intuitively, a set of individual tracks -- there is one track per stellar model, defined uniquely by their initial stellar mass.
 
## Individual Tracks
```@docs
StellarTracks.PARSEC.PARSECTrack
```

## Utilities
In PARSEC the initial helium abundance ``Y`` is scaled with the initial metallicity, allowing for easy conversion between metal mass fraction ``Z`` and the logarithmic metal abundance [M/H].

```@docs
StellarTracks.PARSEC.PARSEC_MH
StellarTracks.PARSEC.PARSEC_Z
```

## PARSEC References
This page cites the following references:

```@bibliography
Pages = ["parsec.md"]
Canonical = false
```
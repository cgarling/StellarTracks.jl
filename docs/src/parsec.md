```@meta
CurrentModule = StellarTracks
ShareDefaultModule = true
```

```@setup
include("plotting.jl")
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

This package handles downloading and pre-processing of the EEP tracks produced by [Rosenfield2016](@citet) (available [here](https://github.com/philrosenfield/padova_tracks)) using [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl). The main access point we provide is [`PARSECLibrary`](@ref StellarTracks.PARSEC.PARSECLibrary), which will load and make available the full library of stellar tracks. The first time you call this method, you will be prompted to download the required data files. The total data volume is ~150 MB after processing. Information on customizing the install location is available [here](https://www.oxinabox.net/DataDeps.jl/stable/z10-for-end-users/). The data can be uninstalled by running `using DataDeps; rm(datadep"PARSECv1.2S"; recursive=true)`. With all the tracks available, we are able to perform operations like interpolating isochrones at any age and metallicity within the PARSEC parameter space.

## Examples
First we load the full PARSEC library, which is downloaded via DataDeps.jl if not already available.
```@example
using StellarTracks.PARSEC
p = PARSECLibrary()
```

Now we use the [`PARSEC.PARSECLibrary`](@ref) to interpolate an isochrone at `log10(age [yr]) = 10.05` and logarithmic metallicity \[M/H\]=-1.234. The isochrone is returned as a `NamedTuple`.
```@example
iso = isochrone(p, 10.05, -1.234)
```

The `NamedTuple` returned by `isochrone` can be converted to table types, like `TypedTables.Table` to simplify further use.
```@example
using TypedTables: Table
Table(iso)
```

The theoretical isochrone is plotted below.

```@example
plot_hr(iso) # hide
```

We can load a grid of bolometric corrections from [BolometricCorrections.jl](@extref BolometricCorrections overview) to add observational magnitudes to the theoretical isochrone. In this example, we use the MIST bolometric correction grid, which offers bolometric corrections for varying metallicities (\[M/H\]) and reddening values (``A_V``).

Because the solar metallicity calibrations of PARSEC and MIST are not exactly the same, the protostellar metal mass fraction ``Z`` that corresponds to a given \[M/H\] is not the same between the two libraries. The `isochrone` interface will convert the given \[M/H\], which is assumed to be the desired metallicity in the *stellar track* library, to its corresponding metal mass fraction, and then convert from the metal mass fraction to the correct \[M/H\] for the assumed chemical model of the bolometric correction grid.

This method returns a `TypedTables.Table` that contains the information from both sources. Here we evaluate an isochrone with `log10(age [yr]) = 10.05`, \[M/H\]=-1.234, and ``A_v=0.02`` mag. 

```@example
using BolometricCorrections.MIST: MISTBCGrid
m = MISTBCGrid("JWST")
iso = isochrone(p, m, 10.05, -1.234, 0.02)
```

All available columns in the isochrone can be obtained with `TypedTables.columnnames`.

```@example
using TypedTables: columnnames
columnnames(iso)
```

A color-magnitude diagram constructed from the isochrone is plotted below.

```@example
plot_cmd(iso) # hide
```

## Chemistry API
We provide the [`StellarTracks.PARSEC.PARSECChemistry`](@ref) type that follows the chemistry API defined in [BolometricCorrections.jl](@extref BolometricCorrections chemistry_api) to access information on the chemical mixture assumed for the PARSEC models.

```@docs
StellarTracks.PARSEC.PARSECChemistry
```

Note that in our conversions between ``Z`` and \[M/H\], remembering that `MH = log10(Z/X) - log10(Z⊙/X⊙)`, we use the *photospheric* solar values for `Z⊙` and `X⊙` (these are `Z_⊙` and `X_⊙ = 1 - Z_⊙ - Y_⊙` in Table 3 of [Bressan2012](@citet)). This reproduces the relation between `Z` and \[M/H\] defined in Table 4 of [Bressan2012](@citet), which is also used in the "CMD" webform provided by the PARSEC team.

## Library API
```@docs
StellarTracks.PARSEC.PARSECLibrary
isochrone(::StellarTracks.PARSEC.PARSECLibrary, ::Number, ::Number)
```

The full library is principally a set of [`PARSECTrackSet`](@ref StellarTracks.PARSEC.PARSECTrackSet)s, with one track set per unique chemical composition. All PARSEC models have scaled-solar chemical compositions, so they vary only in total metallicity (i.e., ``Z``). 

## Track Set API
```@docs
StellarTracks.PARSEC.PARSECTrackSet
```
 
## Individual Tracks API
```@docs
StellarTracks.PARSEC.PARSECTrack
```

## PARSEC References
This page cites the following references:

```@bibliography
Pages = ["parsec.md"]
Canonical = false
```
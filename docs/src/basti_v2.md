```@meta
CurrentModule = StellarTracks
ShareDefaultModule = true
```

```@setup
include("plotting.jl")
```

# [BaSTIv2](@id BaSTIv2)

Here we describe the interface we provide to the updated BaSTI stellar models presented in [Hidalgo2018,Pietrinferni2021,Salaris2022,Pietrinferni2024](@citep). BaSTIv2 specific code is housed in the `BaSTIv2` submodule, which can be accessed as

```@example
using StellarTracks.BaSTIv2 # load all exported methods
using StellarTracks.BaSTIv2: BaSTIv2Library, X, Y, Z # load specific methods
```

## Data Acquisition

The tracks will be downloaded automatically using [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl) the first time you try to access them. The main access point we provide is [`BaSTIv2Library`](@ref StellarTracks.BaSTIv2.BaSTIv2Library), which will load and make available the full library of stellar tracks. The first time you call this method, you will be prompted to download the required data files. The total data volume is **XX** MB. Information on customizing the install location is available [here](https://www.oxinabox.net/DataDeps.jl/stable/z10-for-end-users/). The data can be uninstalled by running `using DataDeps; rm(datadep"BaSTIv2"; recursive=true)`. With all the tracks available, we are able to perform operations like interpolating isochrones at any age and metallicity within the BaSTIv2 parameter space.

## Grid Properties
The BaSTIv2 model grid contains models for the following iron abundances:

```@example
BaSTIv2.feh_grid
```

This is a superset of all available iron abundances; not all combinations of physics and chemistry will have all these iron abundances available. In particular, the α-enhanced models with \[α/Fe\]=0.4 presently have a limited number of iron abundances available.

The BaSTIv2 grid includes models with scaled-solar abundance patterns, α-enhanced models with an average \[α/Fe\]=0.4 (presented in [Pietrinferni2021](@citet)), and α-depleted models with \[α/Fe\]=-0.2 (presented in [Pietrinferni2024](@citet). For these models, the abundances of the α elements O, Ne, Mg, Si, S, Ar, Ca, and Ti have been changed relative to Fe. These α-enhanced models are useful for modeling low-metallicity stars that formed prior to significant iron enrichment from type Ia supernovae. These stars are most common in the Galactic halo and low-mass dwarf galaxies. Note that the conversion between metal mass fraction ``Z`` and logarithmic metal abundance \[M/H\] is the same for the scaled-solar models as for the α-enhanced models, however the iron abundance \[Fe/H\] is not the same as \[M/H\] -- see Table 1 of [Pietrinferni2021](@citet) and [Pietrinferni2024](@citet) for the full elemental abundance tables.

For scaled-solar abundance patterns, models are available with different physics models varying convective overshooting during central H-burning and atomic diffusion. So-called "canonical" models do not including overshooting, while "non-canonical" models do include a convective overshooting treatment (see section 3 of [Pietrinferni2004](@citet)). Inclusion of convective overshooting during central H-burning mainly changes the main sequence turn-off morphology and generally results in better fits to simple stellar populations like globular clusters. For methods in this module that take a `canonical::Bool` argument, a value of `canonical=true` indicates you want to use the "canonical" stellar models, while a value of `canonical=false` means you want to use the "non-canonical" models that include convective core overshooting. Similarly, setting an argument `diffusion::Bool=true` indicates you wish to use the models with atomic diffusion. For scaled-solar abundance patterns, models with different Reimers mass loss parameters `η` are also available.

The α-depleted [Pietrinferni2024](@citep) and α-enhanced [Pietrinferni2021](@citep) models are only available with overshooting and diffusion (`canonical=false`, `diffusion=true`, with `η=0.3`). Additional α-enhanced models with enhanced primordial helium abundances are also available (`yp = 0.247 (fiducial), 0.275, 0.3, 0.32`).

Not all combinations of `canonical`, `diffusion`, α-element abundance, and primordial helium abundance `yp` are valid. The table below summarizes the available combinations of parameters.

| \[α/Fe\] | canonical | diffusion | yp    | η   |
|----------|-----------|-----------|-------|-----|
| -0.2     | false     | true      | 0.247 | 0.3 |
| 0.0      | false     | true      | 0.247 | 0.3 |
| 0.0      | false     | false     | 0.247 | 0.0 |
| 0.0      | false     | false     | 0.247 | 0.3 |
| 0.0      | true      | false     | 0.247 | 0.0 |
| 0.4      | false     | true      | 0.247 | 0.3 |
| 0.4      | false     | true      | 0.275 | 0.3 |
| 0.4      | false     | true      | 0.300 | 0.3 |
| 0.4      | false     | true      | 0.320 | 0.3 |

As setting up these arguments can be arduous, we have attempted to provide useful error messages when an invalid combination of arguments is requested.

## Examples
First we load the non-canonical BaSTIv2 library with diffusion and \[α/Fe\]=0.0, which is downloaded via DataDeps.jl if not already available.
```@example
using StellarTracks.BaSTIv2
p = BaSTIv2Library(0.0, false, true, 0.247, 0.3)
```

Now we use the [`BaSTIv2Library`](@ref) to interpolate an isochrone at `log10(age [yr]) = 10.05` and logarithmic metallicity \[M/H\]=-1.234. The isochrone is returned as a `NamedTuple`.
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

We can load a grid of bolometric corrections from [BolometricCorrections.jl](https://github.com/cgarling/BolometricCorrections.jl) to add observational magnitudes to the theoretical isochrone. In this example, we use the MIST bolometric correction grid, which offers bolometric corrections for varying metallicities (\[M/H\]) and reddening values (``A_V``).

Because the solar metallicity calibrations of BaSTIv2 and MIST are not exactly the same, the protostellar metal mass fraction ``Z`` that corresponds to a given \[M/H\] is not the same between the two libraries. The `isochrone` interface will convert the given \[M/H\], which is assumed to be the desired metallicity in the *stellar track* library, to its corresponding metal mass fraction, and then convert from the metal mass fraction to the correct \[M/H\] for the assumed chemical model of the bolometric correction grid.

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

We provide the [`StellarTracks.BaSTIv2.BaSTIv2Chemistry`](@ref) type that follows the chemistry API defined in [BolometricCorrections.jl](https://github.com/cgarling/BolometricCorrections.jl) to access information on the chemical mixture assumed for the BaSTIv2 models.

```@docs
StellarTracks.BaSTIv2.BaSTIv2Chemistry
```

## Library API
```@docs
StellarTracks.BaSTIv2.BaSTIv2Library
isochrone(::StellarTracks.BaSTIv2.BaSTIv2Library, ::Number, ::Number)
```

The full library is principally a set of [`BaSTIv2TrackSet`](@ref StellarTracks.BaSTIv2.BaSTIv2TrackSet)s, with one track set per unique chemical composition. We do not presently offer interpolation as a function of \[α/Fe\] or between the canonical and non-canonical models, or between the models with and without diffusion, so the individual track sets in the library vary only in total metallicity (i.e., ``Z``).

## Track Set API
```@docs
StellarTracks.BaSTIv2.BaSTIv2TrackSet
```
 
## Individual Tracks API
```@docs
StellarTracks.BaSTIv2.BaSTIv2Track
```

## BaSTIv2 References
This page cites the following references:

```@bibliography
Pages = ["basti_v2.md"]
Canonical = false
```
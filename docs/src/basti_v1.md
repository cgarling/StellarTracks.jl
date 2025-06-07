```@meta
CurrentModule = StellarTracks
ShareDefaultModule = true
```

```@setup
include("plotting.jl")
```

# [BaSTIv1](@id BaSTIv1)

Here we describe the interface we provide to the older BaSTI stellar models presented in [Pietrinferni2004,Pietrinferni2006,Pietrinferni2013](@citet) -- as there are newer BaSTI models (e.g., [Hidalgo2018](@citet)), we describe these older models as `BaSTIv1` to differentiate them. BaSTIv1 specific code is housed in the `BaSTIv1` submodule, which can be accessed as

```@example
using StellarTracks.BaSTIv1 # load all exported methods
using StellarTracks.BaSTIv1: BaSTIv1Library, X, Y, Z # load specific methods
```

## Data Acquisition

The tracks will be downloaded automatically using [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl) the first time you try to access them. The main access point we provide is [`BaSTIv1Library`](@ref StellarTracks.BaSTIv1.BaSTIv1Library), which will load and make available the full library of stellar tracks. The first time you call this method, you will be prompted to download the required data files. The total data volume is ~32 MB. Information on customizing the install location is available [here](https://www.oxinabox.net/DataDeps.jl/stable/z10-for-end-users/). The data can be uninstalled by running `using DataDeps; rm(datadep"BaSTIv1"; recursive=true)`. With all the tracks available, we are able to perform operations like interpolating isochrones at any age and metallicity within the BaSTIv1 parameter space.

## Grid Properties
The BaSTIv1 model grid contains models for the following metal mass fractions:

```@example
BaSTIv1.zgrid
```

which correspond to the following values of \[M/H\]:

```@example
MH.(BaSTIv1Chemistry(), BaSTIv1.zgrid)
```

Note that the BaSTIv1 stellar models only include initial stellar masses from 0.5 to 10 solar masses. They therefore do not include the lower main sequence that can be important for modeling very nearby stellar populations. Nor do they include very high-mass stars (e.g., O-type stars) that can be important when studying populations with high present0-day SFRs.

The BaSTIv1 grid includes models with scaled-solar abundance patterns as well as α-enhanced models with an average \[α/Fe\]=0.4 (presented in [Pietrinferni2006](@citet)). Recall that the α elements are O, Ne, Mg, Si, S, Ar, Ca, and Ti. These α-enhanced models are useful for modeling low-metallicity stars that formed prior to significant iron enrichment from type Ia supernovae. These stars are most common in the Galactic halo and low-mass dwarf galaxies. Note that the conversion between metal mass fraction ``Z`` and logarithmic metal abundance \[M/H\] is the same for the scaled-solar models as for the α-enhanced models, however the iron abundance \[Fe/H\] is not the same as \[M/H\] -- see Table 1 of [Pietrinferni2006](@citet).

The BaSTIv1 models are also differentiated by whether they include convective overshooting during central H-burning. So-called "canonical" models do not including overshooting, while "non-canonical" models do include a convective overshooting treatment (see section 3 of [Pietrinferni2004](@citet)). Inclusion of convective overshooting during central H-burning mainly changes the main sequence turn-off morphology and generally results in better fits to simple stellar populations like globular clusters. For methods in this module that take a `canonical::Bool` argument, a value of `canonical=true` indicates you want to use the "canonical" stellar models, while a value of `canonical=false` means you want to use the "non-canonical" models that include convective core overshooting.

## Examples
First we load the full BaSTIv1 library, which is downloaded via DataDeps.jl if not already available.
```@example
using StellarTracks.BaSTIv1
p = BaSTIv1Library(0.0, false)
```

Now we use the [`BaSTIv1Library`](@ref) to interpolate an isochrone at `log10(age [yr]) = 10.05` and logarithmic metallicity \[M/H\]=-1.234. The isochrone is returned as a `NamedTuple`.
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









## Library API
```@docs
StellarTracks.BaSTIv1.BaSTIv1Library
isochrone(::StellarTracks.BaSTIv1.BaSTIv1Library, ::Number, ::Number)
```

The full library is principally a set of [`BaSTIv1TrackSet`](@ref StellarTracks.BaSTIv1.BaSTIv1TrackSet)s, with one track set per unique chemical composition. We do not presently offer interpolation as a function of \[α/Fe\] or between the canonical and non-canonical models, so the individual track sets in the library vary only in total metallicity (i.e., ``Z``).

## Track Set API
```@docs
StellarTracks.BaSTIv1.BaSTIv1TrackSet
```
 
## Individual Tracks API
```@docs
StellarTracks.BaSTIv1.BaSTIv1Track
```

## BaSTIv1 References
This page cites the following references:

```@bibliography
Pages = ["basti_v1.md"]
Canonical = false
```
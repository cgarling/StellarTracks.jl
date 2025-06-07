```@meta
CurrentModule = StellarTracks
ShareDefaultModule = true
```

```@setup
include("plotting.jl")
```

# [MIST](@id MIST)

Here we describe the interface we provide to the MIST v1.2 library of stellar evolutionary tracks. MIST specific code is housed in the `MIST` submodule, which can be accessed as

```julia
using StellarTracks.MIST # load all exported methods
using StellarTracks.MIST: MISTLibrary, X, Y, Z # load specific methods
```

The main papers describing the MIST family of stellar models are [Dotter2016,Choi2016](@citet). The tracks as provided by the MIST team [here](https://waps.cfa.harvard.edu/MIST/model_grids.html) include the equivalent evolutionary points (EEPs) necessary to support robust isochrone creation and interpolation.

The MIST library has been widely used as it covers the full range of stellar masses and metallicities relevant for most studies of stellar populations. MIST includes stars with initial stellar masses from 0.1 to 300 solar masses and initial metallicities ``-4 \le [\text{M}/\text{H}] \le 0.5``. MIST includes post-MS and post-RGB evolution (when appropriate). MIST also provides rotating (`vvcrit=0.4`) and non-rotating (`vvcrit=0.0`) models.

## Data Acquisition

This package handles downloading and pre-processing of the MIST stellar tracks. The main access point we provide is [`MISTLibrary`](@ref StellarTracks.MIST.MISTLibrary), which will load and make available the full library of stellar tracks. The first time you call this method, you will be prompted to download the required data files. The total data volume that will be downloaded is about 1.3 GB and will total 158 MB after processing. Information on customizing the install location is available [here](https://www.oxinabox.net/DataDeps.jl/stable/z10-for-end-users/). The data can be uninstalled by running `using DataDeps; rm(datadep"MISTv1.2_vvcrit0.0"; recursive=true); rm(datadep"MISTv1.2_vvcrit0.4"; recursive=true)`. With all the tracks available, we are able to perform operations like interpolating isochrones at any age and metallicity within the MIST parameter space.

## Table Details

The user guide for the MIST products is available [here](https://waps.cfa.harvard.edu/MIST/README_overview.pdf). The full MIST tracks contain 77 data columns originating from the MESA output. An description of the columns is available [here](https://waps.cfa.harvard.edu/MIST/README_tables.pdf). **Currently, we process the raw tracks and only save the subset of columns given by `StellarTracks.MIST.select_columns` (see below).** These columns are the ones most commonly needed for computing isochrones and applying bolometric corrections to compare against observed stellar populations. This choice is an optimization for storage space, load time, and development simplicity. If you require access to more columns, please submit an issue on the source repository and we can consider options.
```@example
using StellarTracks.MIST
MIST.select_columns # These columns are saved from raw tracks
```

## Examples
Load the full MIST library of non-rotating models `vvcrit=0`, which is downloaded via DataDeps.jl if not already available. MIST also provides rotating models with `vvcrit=0.4` which can be loaded with `MISTLibrary(0.4)`.
```@example
using StellarTracks.MIST
p = MISTLibrary(0.0)
```

Use the [`MIST.MISTLibrary`](@ref) to interpolate an isochrone at `log10(age [yr]) = 10.05` and metallicity ``[\text{M}/\text{H}] = -1.234``. The isochrone is returned as a `NamedTuple`.
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

We can load a grid of bolometric corrections from [BolometricCorrections.jl](https://github.com/cgarling/BolometricCorrections.jl) to add observational magnitudes to the theoretical isochrone. In this example, we use the MIST bolometric correction grid, which offers bolometric corrections for varying metallicities (\[Fe/H\]) and reddening values (``A_V``). This method returns a `TypedTables.Table` that contains the information from both sources. Here we evaluate an isochrone with `log10(age [yr]) = 10.05`, ``[\text{M}/\text{H}] = -1.234``, and ``A_v=0.02`` mag. 

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
fig,ax1 = plt.subplots() # hide
ax1.plot(iso.F090W .- iso.F150W, iso.F090W) # hide
ax1.set_ylim(reverse(ax1.get_ylim())) # hide
ax1.set_xlim([0.4, 1.62]) # hide
ax1.set_xlabel("F090W - F150W") # hide
ax1.set_ylabel("F090W") # hide
fig # hide
```



## Chemistry API
We re-export the `MIST.MISTChemistry` **add interlink** type defined in BolometricCorrections.jl that can be used to access information on the chemical mixture assumed for the MIST models. All models have scaled-solar chemical compositions.

## Library API

```@docs
StellarTracks.MIST.MISTLibrary
isochrone(::StellarTracks.MIST.MISTLibrary, ::Number, ::Number)
```

## Track Set API
```@docs
StellarTracks.MIST.MISTTrackSet
```
 
## Individual Tracks API
```@docs
StellarTracks.MIST.MISTTrack
```

## PARSEC References
This page cites the following references:

```@bibliography
Pages = ["mist.md"]
Canonical = false
```
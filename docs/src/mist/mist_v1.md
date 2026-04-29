```@meta
CurrentModule = StellarTracks
ShareDefaultModule = true
```

```@setup
include("../plotting.jl")
```

# [MIST v1.2](@id MIST_v1)

Here we describe the interface we provide to the MIST v1.2 library of stellar evolutionary tracks. See the [MIST overview](@ref MIST) page for a comparison with MIST v2.5.

The main papers describing the MIST v1.2 models are [Dotter2016,Choi2016](@citet). The tracks include equivalent evolutionary points (EEPs) necessary to support robust isochrone creation and interpolation.

The MIST v1.2 library covers initial stellar masses from 0.1 to 300 solar masses and initial metallicities ``-4 \le [\text{M}/\text{H}] \le 0.5``. It includes post-MS and post-RGB evolution (when appropriate), and provides both rotating (`vvcrit=0.4`) and non-rotating (`vvcrit=0.0`) models. All models assume scaled-solar chemical compositions using the solar abundances of [Asplund2009](@citet).

## Data Acquisition

This package handles downloading and pre-processing of the MIST stellar tracks. The main access point we provide is [`MISTv1Library`](@ref StellarTracks.MIST.MISTv1Library), which will load and make available the full library of stellar tracks. The first time you call this method, you will be prompted to download the required data files. The total data volume that will be downloaded is about 1.3 GB and will total ~158 MB after processing. Information on customizing the install location is available [here](https://www.oxinabox.net/DataDeps.jl/stable/z10-for-end-users/). The data can be uninstalled by running `using DataDeps; rm(datadep"MISTv1.2_vvcrit0.0"; recursive=true); rm(datadep"MISTv1.2_vvcrit0.4"; recursive=true)`. With all the tracks available, we are able to perform operations like interpolating isochrones at any age and metallicity within the MIST parameter space.

## Table Details

The user guide for the MIST products is available [here](https://mist.science/README_overview.pdf). The full MIST tracks contain 77 data columns originating from the MESA output. A description of the columns is available [here](https://mist.science/README_tables.pdf). **Currently, we process the raw tracks and only save the subset of columns given by `StellarTracks.MIST.select_columns` (see below).** These columns are the ones most commonly needed for computing isochrones and applying bolometric corrections to compare against observed stellar populations. This choice is an optimization for storage space, load time, and development simplicity. If you require access to more columns, please submit an issue on the source repository and we can consider options.

```@example
using StellarTracks.MIST
MIST.select_columns # These columns are saved from raw tracks
```

## Examples
Load the full MIST v1.2 library of non-rotating models (`vvcrit=0`), which is downloaded via DataDeps.jl if not already available. MIST also provides rotating models with `vvcrit=0.4` which can be loaded with `MISTv1Library(0.4)`.

```@example
using StellarTracks.MIST
p = MISTv1Library(0.0)
```

Use the [`MISTv1Library`](@ref StellarTracks.MIST.MISTv1Library) to interpolate an isochrone at `log10(age [yr]) = 10.05` and metallicity ``[\text{M}/\text{H}] = -1.234``. The isochrone is returned as a `NamedTuple`.

```@example
iso = isochrone(p, 10.05, -1.234)
```

The `NamedTuple` returned by `isochrone` can be converted to table types, like `TypedTables.Table`, to simplify further use.

```@example
using TypedTables: Table
Table(iso)
```

The theoretical isochrone is plotted below.

```@example
plot_hr(iso) # hide
```

We can load a grid of bolometric corrections from [BolometricCorrections.jl](@extref BolometricCorrections overview) to add observational magnitudes to the theoretical isochrone. BolometricCorrections.jl provides two versions of the MIST BC grid: v1.2 (metallicity and reddening only) and v2.5 (adds \[α/Fe\]). Here we use the **MIST v1.2** grid to evaluate an isochrone with `log10(age [yr]) = 10.05`, ``[\text{M}/\text{H}] = -1.234``, and ``A_v=0.02`` mag.

```@example
using BolometricCorrections.MIST: MISTv1BCGrid
m1 = MISTv1BCGrid("JWST")
iso_v1 = isochrone(p, m1, 10.05, -1.234, 0.02)
```

All available columns in the isochrone can be obtained with `TypedTables.columnnames`.

```@example
using TypedTables: columnnames
columnnames(iso_v1)
```

We can also use the **MIST v2.5** BC grid. The `isochrone` method accepts a `MISTv2BCGrid` and converts metallicities between the chemical abundance scales of the stellar tracks and bolometric corrections automatically. The \[α/Fe\] value is inferred from the track library automatically.

```@example
using BolometricCorrections.MIST: MISTv2BCGrid
m2 = MISTv2BCGrid("JWST")
iso_v2 = isochrone(p, m2, 10.05, -1.234, 0.02) # [M/H], Av
```

A color-magnitude diagram comparing the two BC grids applied to an isochrone made from MIST v1.2 tracks is plotted below. The isochrone made with the v1.2 BCs is shown as a solid line and the isochrone made with the v2.5 BCs is shown as a dashed line.

```@example
colors_v1 = getproperty(iso_v1, :F090W) .- getproperty(iso_v1, :F150W) # hide
mags_v1 = getproperty(iso_v1, :F090W) # hide
colors_v2 = getproperty(iso_v2, :NIRCAM_F090W) .- getproperty(iso_v2, :NIRCAM_F150W) # hide
mags_v2 = getproperty(iso_v2, :NIRCAM_F090W) # hide
fig = Figure(size=(500, 500)) # hide
ax = Axis(fig[1,1], # hide
            xlabel="F090W - F150W", # hide
            ylabel="F090W", # hide
            yreversed=true, # hide
            limits=(0.4, 1.62, nothing, nothing)) # hide
lines!(ax, colors_v1, mags_v1; label="v1") # hide
lines!(ax, colors_v2, mags_v2; linestyle=:dash, label="v2") # hide
axislegend(ax, position=:rc) # hide
fig # hide
```

## Chemistry API
We re-export the [`BolometricCorrections.MIST.MISTv1Chemistry`](@extref) type defined in BolometricCorrections.jl that can be used to access information on the chemical mixture assumed for the MIST v1.2 models. All models have scaled-solar chemical compositions.

## Library API

```@docs
StellarTracks.MIST.MISTv1Library
isochrone(::StellarTracks.MIST.MISTv1Library, ::Number, ::Number)
```

## Track Set API
```@docs
StellarTracks.MIST.MISTv1TrackSet
```

## Individual Tracks API
```@docs
StellarTracks.MIST.MISTv1Track
```

## MIST v1.2 References
This page cites the following references:

```@bibliography
Pages = ["mist_v1.md"]
Canonical = false
```

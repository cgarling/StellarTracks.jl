```@meta
CurrentModule = StellarTracks
ShareDefaultModule = true
```

```@setup
import PyPlot as plt
plt.ioff()
ENV["MPLBACKEND"] = "agg"
import PyPlot: @L_str # For LatexStrings
plt.rc("text", usetex=true)
plt.rc("font", family="serif", serif=["Computer Modern"], size=16)
# This gets close but not quite
# plt.matplotlib.rcParams["axes.formatter.use_mathtext"] = true
# plt.rc("font", family="serif", serif=["cmr10"], size=14)
plt.rc("figure", figsize=(5,5))
plt.rc("patch", linewidth=1, edgecolor="k", force_edgecolor=true)
# https://matplotlib.org/stable/gallery/images_contours_and_fields/interpolation_methods.html
plt.rc("image", interpolation="none")
```

# [MIST](@id MIST)

Here we describe the interface we provide to the MIST v1.2 library of stellar evolutionary tracks. MIST specific code is housed in the `MIST` submodule, which can be accessed as

```julia
using StellarTracks.MIST # load all exported methods
using StellarTracks.MIST: MISTLibrary, X, Y, Z # load specific methods
```

The main papers describing the MIST family of stellar models are [Dotter2016,Choi2016](@citet). The tracks as provided by the MIST team [here](https://waps.cfa.harvard.edu/MIST/model_grids.html) include the equivalent evolutionary points (EEPs) necessary to support robust isochrone creation and interpolation.

## Data Acquisition

This package handles downloading and pre-processing of the MIST stellar tracks. The main access point we provide is [`MISTLibrary`](@ref StellarTracks.MIST.MISTLibrary), which will load and make available the full library of stellar tracks. The first time you call this method, you will be prompted to download the required data files. The total data volume that will be downloaded is about 1.3 GB and will total **XX GB** after additional processing. Information on customizing the install location is available [here](https://www.oxinabox.net/DataDeps.jl/stable/z10-for-end-users/). The data can be uninstalled by running `using DataDeps; rm(datadep"PARSECv1.2S"; recursive=true)`. With all the tracks available, we are able to perform operations like interpolating isochrones at any age and metallicity within the MIST parameter space.

## Examples
Load the full MIST library, which is downloaded via DataDeps.jl if not already available.
```@example
using StellarTracks.MIST
p = MISTLibrary()
```

Use the [`MIST.MISTLibrary`](@ref) to interpolate an isochrone at `log10(age [yr]) = 10.05` and metal mass fraction ``Z=0.001654``. The isochrone is returned as a `NamedTuple`.
```@example
iso = isochrone(p, 10.05, 0.001654)
```

## Table Details
The user guide for the MIST products is available [here](https://waps.cfa.harvard.edu/MIST/README_overview.pdf). The full MIST tracks contain 77 data columns originating from the MESA output. An description of the columns is available [here](https://waps.cfa.harvard.edu/MIST/README_tables.pdf). 
```@meta
CurrentModule = StellarTracks
ShareDefaultModule = true
```

```@setup mistv2
include("../plotting.jl")
```

# [MIST v2.5](@id MIST_v2)

Here we describe the interface we provide to the MIST v2.5 library of stellar evolutionary tracks. See the [MIST overview](@ref MIST) page for a comparison with MIST v1.2.

The MIST v2.5 tracks are described in [Dotter2026](@citet) (α-enhanced models) and [Bauer2026](@citet) (white dwarf cooling sequences). The key advances relative to [MIST v1.2](@ref MIST_v1) are:

- **\[α/Fe\] as a free parameter**: five \[α/Fe\] values spanning −0.2 to +0.6 are provided, making MIST v2.5 particularly well-suited for studying metal-poor stellar populations where α-enhancement is common.
- **Revised solar abundances**: the v2.5 models adopt solar abundances from [Grevesse1998](@citet) rather than [Asplund2009](@citet) used in v1.2. The `MISTv2Chemistry` type (provided by BolometricCorrections.jl) encodes this mixture.
- **Denser metallicity grid**: 17 \[Fe/H\] values (vs. 15 in v1.2), adding −2.75 and −2.25 dex.
- **Altered initial mass grid**: 156 mass points (vs. 196 in v1.2) but with denser sampling in the range 0.7–1.2 M☉.

## Data Acquisition

The MIST v2.5 tracks are organized by rotation parameter (`vvcrit`) and α-element enhancement (`afe`). Each (`vvcrit`, `afe`) combination is a separate DataDep. The main access point is [`MISTv2Library`](@ref StellarTracks.MIST.MISTv2Library), which takes both parameters as arguments:

```julia
lib = MISTv2Library(0.0, 0.0)   # non-rotating, solar-scaled [α/Fe]=0
lib = MISTv2Library(0.4, 0.4)   # rotating, α-enhanced [α/Fe]=+0.4
```

The available values are:
- `vvcrit`: `StellarTracks.MIST.vvcrit_grid_v2` = {0.0, 0.4}
- `afe`: `StellarTracks.MIST.afe_grid_v2` = {−0.2, 0.0, 0.2, 0.4, 0.6}

The first time you call `MISTv2Library` for a given (`vvcrit`, `afe`) combination, you will be prompted to download the required data files. Each combination requires downloading 16–17 `.txz` archives (one per \[Fe/H\] value; the \[α/Fe\] = +0.6 grid has 16 values because \[Fe/H\] = +0.5 is not available). The total data volume per combination is approximately 1.3–1.5 GB before processing and ~160 MB after. Information on customizing the install location is available [here](https://www.oxinabox.net/DataDeps.jl/stable/z10-for-end-users/).

!!! note "Data layout"
    Each downloaded DataDep is named `MISTv2.5_vvcrit<V>_afe_<A>`, e.g. `MISTv2.5_vvcrit0.0_afe_p0` for `vvcrit=0.0, afe=0.0`. To remove a specific combination, run e.g. `using DataDeps; rm(datadep"MISTv2.5_vvcrit0.0_afe_p0"; recursive=true)`.

## Table Details

As with v1.2, only the subset of columns given by `StellarTracks.MIST.select_columns` is saved after processing. The MIST v2.5 tracks use the same column format as v1.2.

```@example mistv2
using StellarTracks.MIST
MIST.select_columns
```

Note that the MIST v2.5 grid has some individual mass tracks missing for certain (`feh`, `vvcrit`, `afe`) combinations. `MISTv2TrackSet` silently skips missing files rather than erroring, so the set of available masses may differ across the parameter grid.

## Examples

Load the MIST v2.5 library of non-rotating, solar-scaled models. The `vvcrit` and `afe` parameters must be exact values from the grid.

```@example mistv2
using StellarTracks.MIST
lib = MISTv2Library(0.0, 0.0)
```

Interpolate an isochrone at `log10(age [yr]) = 10.05` and ``[\text{M}/\text{H}] = -1.234``. The isochrone is returned as a `NamedTuple`.

```@example mistv2
iso = isochrone(lib, 10.05, -1.234)
```

The result can be converted to a table:

```@example mistv2
using TypedTables: Table
Table(iso)
```

Plot the HR diagram:

```@example mistv2
plot_hr(iso) # hide
```

Load a bolometric correction grid from [BolometricCorrections.jl](@extref BolometricCorrections overview) to place the isochrone into an observational photometric system. Because the MIST v2.5 stellar track and BC grids share the same `MISTv2Chemistry` chemical mixture, metallicities are handled consistently without the need for any conversion.

```@example mistv2
using BolometricCorrections.MIST: MISTv2BCGrid
bc = MISTv2BCGrid("JWST")
iso_phot = isochrone(lib, bc, 10.05, -1.234, 0.02)  # logAge, [M/H], Av
```

```@example mistv2
using TypedTables: columnnames
columnnames(iso_phot)
```

To take advantage of the \[α/Fe\] degree of freedom, load a library with a non-zero `afe` value. Here we compare solar-scaled and α-enhanced isochrones at the same age and \[M/H\]. Note that both `isochrone` calls receive the same value of \[M/H\], not \[Fe/H\]; at the same \[M/H\] the α-enhanced population has a lower \[Fe/H\] but the same total metal mass fraction `Z`.

```@example mistv2
lib_afe = MISTv2Library(0.0, 0.4)
iso_afe = isochrone(lib_afe, 10.05, -1.234)
fig = Figure(size=(500, 500)) # hide
ax = Axis(fig[1,1], xlabel="logTe", ylabel="Mbol", # hide
          xreversed=true, yreversed=true, # hide
          limits=(3.5, 3.85, nothing, nothing)) # hide
lines!(ax, iso.logTe, iso.Mbol; label="[α/Fe]=0.0") # hide
lines!(ax, iso_afe.logTe, iso_afe.Mbol; linestyle=:dash, label="[α/Fe]=+0.4") # hide
axislegend(ax, position=:rb) # hide
fig # hide
```

## Chemistry API

We re-export the [`BolometricCorrections.MIST.MISTv2Chemistry`](@extref) type defined in BolometricCorrections.jl that encodes the solar chemical mixture assumed for the MIST v2.5 models ([Grevesse1998](@citet)). Because \[α/Fe\] is a free parameter in v2.5, the data grid is indexed by \[Fe/H\] rather than \[M/H\]. However, the `isochrone` interface always accepts \[M/H\] (the total logarithmic metal abundance), and `MISTv2Chemistry` handles the conversion between \[M/H\] and \[Fe/H\] via the [Salaris1993](@citet) relation. When a `MISTv2Library` is used with a `MISTv2BCGrid`, the `afe` value is inferred from the track library automatically so the correct bolometric corrections are selected.

## Library API

```@docs
StellarTracks.MIST.MISTv2Library
isochrone(::StellarTracks.MIST.MISTv2Library, ::Number, ::Number)
```

## Track Set API

```@docs
StellarTracks.MIST.MISTv2TrackSet
```

## Individual Tracks API

```@docs
StellarTracks.MIST.MISTv2Track
```

## MIST v2.5 References
This page cites the following references:

```@bibliography
Pages = ["mist_v2.md"]
Canonical = false
```

```@meta
CurrentModule = StellarTracks
ShareDefaultModule = true
```

```@setup overview_cmd
include("../plotting.jl")
```

# [MIST](@id MIST)

The [Mesa Isochrones and Stellar Tracks](https://mist.science) (MIST; [Dotter2016,Choi2016,Dotter2026,Bauer2026](@citet)) library provides stellar evolutionary tracks for a wide range of initial stellar masses and metallicities, computed with the [MESA](https://docs.mesastar.org) stellar evolution code. The tracks use the equivalent evolutionary point (EEP) framework [Dotter2016](@citet) to facilitate robust isochrone construction and mass interpolation.

MIST-specific code lives in the `MIST` submodule, which can be accessed as

```julia
using StellarTracks.MIST # load all exported methods
using StellarTracks.MIST: MISTv1Library, X, Y, Z # load specific methods
```

Two versions of the MIST stellar track libraries are currently supported:

- **[MIST v1.2](@ref MIST_v1)** — the original release with scaled-solar chemical compositions, covering initial stellar masses from 0.1 to 300 M☉ and metallicities ``-4 \le [\text{M}/\text{H}] \le 0.5`` assuming the [Asplund2009](@citet) solar chemical abundances.
- **[MIST v2.5](@ref MIST_v2)** — updated release adding \[α/Fe\] as a free parameter, revised solar chemical abundances [Grevesse1998](@cite), a denser metallicity grid, and an expanded initial stellar mass grid [Dotter2026,Bauer2026](@cite).

## Common features

Both versions share the following characteristics:

- Stellar evolutionary tracks computed using MESA, spanning from the pre-main sequence through (when applicable) the white dwarf cooling sequence.
- Both rotating (`vvcrit=0.4`) and non-rotating (`vvcrit=0.0`) models.
- Tracks with equivalent evolutionary points (EEPs) for robust isochrone construction.
- The same subset of data columns extracted from the full tracks.
- Integration with [BolometricCorrections.jl](@extref BolometricCorrections overview) to generate photometric isochrones in any supported bandpass.

## Version comparison

| Feature | v1.2 | v2.5 |
|---------|------|------|
| Primary reference(s) | [Dotter2016,Choi2016](@citet) | [Dotter2026,Bauer2026](@citet) |
| Solar abundances | [Asplund2009](@citet) | [Grevesse1998](@citet) |
| α-element enhancement | fixed (scaled-solar) | \[α/Fe\] ∈ {−0.2, 0, 0.2, 0.4, 0.6} |
| \[Fe/H\] range | −4.0 to +0.5 (15 values) | −4.0 to +0.5 (17 values) |
| Initial mass range | 0.1 to 300 M☉ | 0.1 to 300 M☉ |
| `vvcrit` values | {0.0, 0.4} | {0.0, 0.4} |
| Chemistry type | [`MISTv1Chemistry`](@extref BolometricCorrections.MIST.MISTv1Chemistry) | [`MISTv2Chemistry`](@extref BolometricCorrections.MIST.MISTv2Chemistry) |
| Track type | [`MISTv1Track`](@ref StellarTracks.MIST.MISTv1Track) | [`MISTv2Track`](@ref StellarTracks.MIST.MISTv2Track) |
| Track set type | [`MISTv1TrackSet`](@ref StellarTracks.MIST.MISTv1TrackSet) | [`MISTv2TrackSet`](@ref StellarTracks.MIST.MISTv2TrackSet) |
| Library type | [`MISTv1Library`](@ref StellarTracks.MIST.MISTv1Library) | [`MISTv2Library`](@ref StellarTracks.MIST.MISTv2Library) |

The color-magnitude diagram below illustrates the effect of switching from MIST v1.2 tracks + BCs to MIST v2.5 tracks + BCs for the same age (`log10(age [yr]) = 10.05`), metallicity (``[\text{Fe}/\text{H}] = -1.234``), and reddening (``A_V = 0.02`` mag). Both use `vvcrit=0` and `[α/Fe]=0`. The v1.2 isochrone is shown as a solid line; the v2.5 isochrone as a dashed line.

```@example overview_cmd
using StellarTracks.MIST # hide
using BolometricCorrections.MIST: MISTv1BCGrid, MISTv2BCGrid # hide
p_v1 = MISTv1Library(0.0) # hide
p_v2 = MISTv2Library(0.0, 0.0) # hide
bc_v1 = MISTv1BCGrid("JWST") # hide
bc_v2 = MISTv2BCGrid("JWST") # hide
iso_v1 = isochrone(p_v1, bc_v1, 10.05, -1.234, 0.02) # hide
iso_v2 = isochrone(p_v2, bc_v2, 10.05, -1.234, 0.0, 0.02) # hide
colors_v1 = getproperty(iso_v1, :F090W) .- getproperty(iso_v1, :F150W) # hide
mags_v1   = getproperty(iso_v1, :F090W) # hide
colors_v2 = getproperty(iso_v2, :NIRCAM_F090W) .- getproperty(iso_v2, :NIRCAM_F150W) # hide
mags_v2   = getproperty(iso_v2, :NIRCAM_F090W) # hide
fig = Figure(size=(500, 500)) # hide
ax = Axis(fig[1,1], # hide
          xlabel="F090W − F150W", # hide
          ylabel="F090W", # hide
          yreversed=true, # hide
          limits=(0.4, 1.62, nothing, 13.0)) # hide
lines!(ax, colors_v1, mags_v1; label="v1.2 tracks and BCs") # hide
lines!(ax, colors_v2, mags_v2; linestyle=:dash, label="v2.5 tracks and BCs") # hide
axislegend(ax, position=:rc) # hide
fig # hide
```

Note that for ease of transition from earlier versions of StellarTracks.jl that only supported MIST v1.2, the following deprecated aliases are available:

- `MISTTrack` aliases for [`MISTv1Track`](@ref StellarTracks.MIST.MISTv1Track)
- `MISTTrackSet` aliases for [`MISTv1TrackSet`](@ref StellarTracks.MIST.MISTv1TrackSet)
- `MISTLibrary` aliases for [`MISTv1Library`](@ref StellarTracks.MIST.MISTv1Library)

## MIST References
This page cites the following references:

```@bibliography
Pages = ["mist_overview.md"]
Canonical = false
```

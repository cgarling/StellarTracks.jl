# StellarTracks.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cgarling.github.io/StellarTracks.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cgarling.github.io/StellarTracks.jl/dev/)
[![Build Status](https://github.com/cgarling/StellarTracks.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cgarling/StellarTracks.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cgarling/StellarTracks.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cgarling/StellarTracks.jl)

This package provides access to and interpolation of pre-computed libraries of stellar tracks. Secondary operations like isochrone creation are also supported. Currently included libraries are

 - [PARSECv1.2S](https://stev.oapd.inaf.it/PARSEC/papers.html)
 - [MIST](https://waps.cfa.harvard.edu/MIST/index.html)
 - [BaSTIv1](http://basti-iac.oa-abruzzo.inaf.it/basti-old/) (older models, circa 2013)
 - [BaSTIv2](http://basti-iac.oa-abruzzo.inaf.it/index.html) (updated BaSTI tracks)

These tracks contain only quantities from stellar interior modeling (e.g., bolometric luminosities and effective temperatures) and must be combined with bolometric corrections to make predictions for observations in different filters/bandpasses. This package integrates with [BolometricCorrections.jl](https://github.com/cgarling/BolometricCorrections.jl) to support isochrone interpolation with the addition of bolometric corrections. See the examples below for a quick overview and our documentation linked in the badges above for more details.

## Example Usage
The most common use cases for this package will involve interacting with subtypes of `AbstractTrackLibrary` -- these subtypes (`PARSECLibrary, MISTLibrary, BaSTIv1Library, BaSTIv2Library`) load all relevant stellar tracks for a particular set of model input parameters -- for example, MIST offers both rotating and non-rotating stellar models, and users can choose between these sets when constructing a `MISTLibrary`. These types provide the ability to interpolate tracks to different stellar initial masses and metallicities and interpolate isochrones at different ages and metallicities.

```julia
julia> using StellarTracks

# Load non-rotating models (vvcrit=0)
julia> tracklib = MISTLibrary(0)

# Interpolate a new track at [M/H] = -2.09, M_ini = 1.09 M⊙
julia> newtrack = tracklib(-2.09, 1.09)
InterpolatedTrack with M_ini=1.09, MH=-2.09, Z=0.00012111671593180112, Y=0.2491818874546179, X=0.7506969958294503.

# Interpolate newtrack at log(age [yr]) = 9
julia> newtrack(9)
(log_L = 0.4956590459743111, log_Teff = 3.898620504847251, log_g = 4.528985420811569, log_surf_cell_z = -4.641158732471512)

# Interpolate isochrone at [M/H] = -2.09, log(age [yr]) = 9
julia> iso = isochrone(tracklib, 9, -2.09);

# `isochrone` returns a `NamedTuple` of vectors that can be parsed
# as a table, for example `TypedTables.Table` or `DataFrames.DataFrame`.
julia> using TypedTables: Table

julia> Table(iso)
Table with 7 columns and 609 rows:
      eep  m_ini     logTe    Mbol     logg     logL      log_surf_cell_z
    ┌────────────────────────────────────────────────────────────────────
 1  │ 200  0.100631  3.54677  11.6724  5.35538  -2.77296  -3.9351
 2  │ 201  0.101926  3.54824  11.6322  5.35036  -2.7569   -3.9351
 3  │ 202  0.103545  3.55008  11.5823  5.34418  -2.73691  -3.9351
 ...
```

This package integrates with [BolometricCorrections.jl](https://github.com/cgarling/BolometricCorrections.jl) to support isochrone interpolation with the addition of bolometric corrections.

```julia
julia> using BolometricCorrections

# Load MIST BCs for HST ACS/WFC
julia> bcgrid = MISTBCGrid("hst_acs_wfc")
MIST bolometric correction grid for photometric system MIST_HST_ACS_WFC

# MIST BCs support different reddening (Av)
# Interpolate isochrone at [M/H] = -2.09, log(age [yr]) = 9, Av = 0.13 mag
julia> isochrone(tracklib, bcgrid, 9, -2.09, 0.13)
Table with 20 columns and 609 rows:
      eep  m_ini     logTe    Mbol     logg     logL      log_surf_cell_z  ACS_WFC_F435W  ACS_WFC_F475W  ACS_WFC_F502N  ACS_WFC_F550M  ⋯
    ┌───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 1  │ 200  0.100631  3.54677  11.6724  5.35538  -2.77296  -3.9351          14.7118        14.0111        14.0181        12.7635        ⋯
 2  │ 201  0.101926  3.54824  11.6322  5.35036  -2.7569   -3.9351          14.6536        13.9571        13.9586        12.7166        ⋯
 3  │ 202  0.103545  3.55008  11.5823  5.34418  -2.73691  -3.9351          14.5806        13.8895        13.8844        12.6579        ⋯
```
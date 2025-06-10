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

These tracks contain only quantities from stellar interior modeling (e.g., bolometric luminosities and effective temperatures) and must be combined with bolometric corrections to make predictions for observations in different filters/bandpasses. This package integrates with [BolometricCorrections.jl](https://github.com/cgarling/BolometricCorrections.jl) to support isochrone interpolation with the addition of bolometric corrections. See our documentation linked in the badges above for more details.
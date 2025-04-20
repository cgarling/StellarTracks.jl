# StellarTracks.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cgarling.github.io/StellarTracks.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cgarling.github.io/StellarTracks.jl/dev/)
[![Build Status](https://github.com/cgarling/StellarTracks.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cgarling/StellarTracks.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cgarling/StellarTracks.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cgarling/StellarTracks.jl)

Provides access to and interpolation of pre-computed libraries of stellar tracks. Currently supported libraries are

 - [PARSECv1.2S](https://stev.oapd.inaf.it/PARSEC/papers.html)

These tracks contain only quantities from stellar interior modeling (e.g., bolometric luminosities and effective temperatures) and must be combined with bolometric corrections to make predictions for observations in different filters/bandpasses. This package integrates with [BolometricCorrections.jl](https://github.com/cgarling/BolometricCorrections.jl) to support these use cases. 
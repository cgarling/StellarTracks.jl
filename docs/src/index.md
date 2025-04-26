```@meta
CurrentModule = StellarTracks
```

# [Overview](@id overview)

The purpose of this package is to provide access to pre-computed libraries of stellar evolutionary tracks. These libraries provide predictions for stellar interior and photosphere properties for stars of different initial masses and chemical compositions. In their standard form, each track describes the evolution of a single model star as it evolves through time. Tracks are identified uniquely by their initial stellar masses (``M_\text{ini}``) and chemical compositions. Most libraries will provide, at minimum, a set of stellar models with different total metallicity (typically quantified by the metal mass fraction ``Z``). Some libraries will also vary other parameters like α-element abundance.

Isochrones representing stellar populations of uniform age and chemical composition but varying initial stellar masses can be formed by interpolating the tracks with the correct chemical compositions to the requested age. This form of isochrone creation can struggle when the spacing in initial mass between the tracks in a set is poor and generally does a poor job of capturing isochrone features that only manifest over a small range of initial stellar masses for a given age. Isochrone creation can be improved by identifying equivalent evolutionary points (EEPs) in these tracks [Dotter2016](@citep); example EEPs are the main sequence turn-off and the tip of the red giant branch. For each EEP, a relation between the stellar initial mass and the age of the star at the EEP can be constructed so that the initial mass of the star as a function of age can be inferred for each EEP. Properties like bolometric luminosity can then be interpolated and evaluated as a function of initial mass for this EEP. When there are many well-spaced EEPs, this method of isochrone creation outperforms the easier approach. *We focus on tracks with EEPs to support this use case*, as **isochrone generation is one of the primary features of this package.**

## Supported Libraries
We currently support the following stellar track libraries:
 - [PARSECv1.2S](@ref PARSEC)

## Chemistry
We include information on the chemical mixtures assumed in each supported library above. We use the interface defined in [BolometricCorrections.jl](https://cgarling.github.io/BolometricCorrections.jl/stable/) to provide this information.

## Bolometric Corrections
The core function of this package is to provide an interface to libraries of stellar evolutionary tracks. One of the most common uses of these tracks is to make predictions for stellar observations, which requires placing the theoretical tracks (which model only stellar interior evolution) into the observed space by applying model stellar atmospheres.

For integrated stellar photometry, this is typically done by applying [bolometric corrections](https://en.wikipedia.org/wiki/Bolometric_correction) that integrate a stellar atmosphere over an observational bandpass to predict the stellar model's magnitude in that bandpass. An interface to these bolometric correction grids is provided by [BolometricCorrections.jl](https://github.com/cgarling/BolometricCorrections.jl). We take advantage of this interface to support the generation of isochrones in observational bandpasses.

## API
Our full API documentation is available [here](@ref api).

## References
This page cites the following references:

```@bibliography
Pages = ["index.md"]
Canonical = false
```
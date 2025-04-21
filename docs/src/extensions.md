# [Extensions](@id extensions)

The core function of this package is to provide an interface to libraries of stellar evolutionary tracks. One of the most common uses of these tracks is to make predictions for stellar observations, which requires placing the theoretical tracks (which model only stellar interior evolution) into the observed space by applying model stellar atmospheres.

For integrated stellar photometry, this is typically done by applying [bolometric corrections](https://en.wikipedia.org/wiki/Bolometric_correction) that integrate a stellar atmosphere over an observational bandpass to predict the stellar model's magnitude in that bandpass. An interface to these bolometric correction grids is provided by [BolometricCorrections.jl](https://github.com/cgarling/BolometricCorrections.jl) and we provide additional functionality in this package if BolometricCorrections.jl is also loaded through a package extension.

## [BolometricCorrections.jl](@id bcs)

One of the most common products observers require are isochrones to predict where stars of different initial masses but identical chemical compositions should lie in an observed color-magnitude diagram. We already provide the [`isochrone`](@ref) method to calculate isochrones from the stellar evolutionary tracks, and we offer additional call signatures when a grid of bolometric corrections is available to predict magnitudes in the observational space. Examples of this use case are given on the [PARSEC](@ref PARSEC) page.
# Derivation of the solar alpha-element mass fraction for BaSTIv2Chemistry.
#
# The updated BaSTI stellar evolution models (Hidalgo et al. 2018, Pietrinferni et al.
# 2021, Salaris et al. 2022, Pietrinferni et al. 2024) adopt the solar photospheric
# abundances from:
#   Caffau, E. et al. (2011), Solar Physics, 268, 255–269
#   DOI: 10.1007/s11207-010-9541-4
#   ADS: https://ui.adsabs.harvard.edu/abs/2011SoPh..268..255C
#   (hereafter Caffau2011)
#
# For the alpha-enhanced BaSTIv2 models, the alpha elements that are uniformly modified
# with respect to the Fe abundance are: O, Ne, Mg, Si, S, Ca, Ti.
# Notably, Ar is NOT among the modified alpha elements in BaSTIv2 (unlike BaSTIv1 and
# MIST/ATLAS9/PARSEC), and is therefore excluded from the alpha element set when
# calculating the Salaris et al. (1993) correction factor.
#
# The alpha-element mass fraction f_α is defined as the fraction of total metal (Z) mass
# contributed by alpha elements: O, Ne, Mg, Si, S, Ca, Ti.
# It enters the Salaris et al. (1993) conversion between [M/H] and [Fe/H]:
#   [M/H] = [Fe/H] + log10(f_α · 10^[α/Fe] + (1 - f_α))
#
# For each element i, the mass contribution is proportional to A_i · 10^(ε_i - 12),
# where A_i is the atomic mass number and ε_i = log10(N_i / N_H) + 12 is the standard
# logarithmic abundance on the scale where log(N_H) + 12 = 12.
#
# Representative metals from Table B.1 of Caffau2011 are listed below, sufficient to
# evaluate f_α.  Omitted trace elements contribute negligibly to Z and do not affect
# f_α at the precision retained here.

using Test: @test
using StellarTracks.BaSTIv2: BaSTIv2Chemistry, alpha_mass_fraction

# ────────────────────────────────────────────────────────────────────────────
# Caffau et al. (2011) Table B.1 photospheric log-epsilon abundances: ε_i = log(N_i/N_H)+12
# Columns: element => (atomic mass number A, log_epsilon)
# Ar is listed for completeness but is NOT an alpha element in BaSTIv2.
# ────────────────────────────────────────────────────────────────────────────
const _caffau2011 = (
    H  = ( 1, 12.00),
    He = ( 4, 10.93),
    C  = (12,  8.50),
    N  = (14,  7.86),
    O  = (16,  8.76),   # alpha
    Ne = (20,  7.84),   # alpha
    Na = (23,  6.29),
    Mg = (24,  7.53),   # alpha
    Al = (27,  6.43),
    Si = (28,  7.51),   # alpha
    P  = (31,  5.41),
    S  = (32,  7.12),   # alpha
    Ar = (36,  6.40),   # NOT alpha for BaSTIv2 (not modified in alpha-enhanced models)
    Ca = (40,  6.31),   # alpha
    Ti = (48,  4.90),   # alpha
    Fe = (56,  7.52),
    Ni = (58,  6.23),
)

const _alpha_elements_caffau2011 = (:O, :Ne, :Mg, :Si, :S, :Ca, :Ti)  # Ar excluded

# Mass contribution of each element: A_i · 10^(ε_i − 12)
_mass(A, ε) = A * exp10(ε - 12.0)

_contributions_caffau2011 = NamedTuple{keys(_caffau2011)}(
    _mass(A, ε) for (A, ε) in values(_caffau2011)
)

Z_total_caffau2011 = sum(v for (el, v) in pairs(_contributions_caffau2011) if el ∉ (:H, :He))
Z_alpha_caffau2011 = sum(v for (el, v) in pairs(_contributions_caffau2011)
                         if el ∈ _alpha_elements_caffau2011)

# f_α = fraction of total metal mass contributed by alpha elements
f_alpha_caffau2011 = Z_alpha_caffau2011 / Z_total_caffau2011

# The stored value alpha_mass_fraction(BaSTIv2Chemistry()) = 0.6475 is f_alpha rounded
# to four decimal places.  We confirm the calculation is consistent with the stored value.
# BaSTIv2Chemistry is parametric (takes α_fe and yp), but alpha_mass_fraction is
# independent of those parameters, so we use representative solar-scaled values.
@test round(f_alpha_caffau2011; digits=4) == alpha_mass_fraction(BaSTIv2Chemistry(0.0, 0.247))

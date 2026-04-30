# Derivation of the solar alpha-element mass fraction for BaSTIv1Chemistry.
#
# The older BaSTI stellar evolution models (Pietrinferni et al. 2004, 2006, 2013) adopt
# the heavy-element distribution of:
#   Grevesse, N. & Noels, A. (1993), in "Origin and Evolution of the Elements",
#   eds. Prantzos, Vangioni-Flam & Cassé, Cambridge University Press, p. 15
#   (hereafter GN93)
#
# The docstring of BaSTIv1Chemistry notes that the GN93 distribution is "minimally
# different" from Grevesse & Sauval (1998) [GS98].  The most significant difference
# is the iron abundance: GN93 gives log ε(Fe) = 7.67 vs. GS98's 7.50.
#
# The alpha-element mass fraction f_α is defined as the fraction of total metal (Z) mass
# contributed by alpha elements: O, Ne, Mg, Si, S, Ca, Ti.
# Note: Ar is NOT counted as an enhanced alpha element in BaSTIv1, consistent with
# BaSTIv2.  Including Ar would give f_α ≈ 0.6682 (rounds to 0.6682), which does not
# match the stored value.
#
# f_α enters the Salaris et al. (1993) conversion between [M/H] and [Fe/H]:
#   [M/H] = [Fe/H] + log10(f_α · 10^[α/Fe] + (1 - f_α))
#
# For each element i, the mass contribution is proportional to A_i · 10^(ε_i - 12),
# where A_i is the atomic mass number and ε_i = log10(N_i / N_H) + 12 is the standard
# logarithmic abundance on the scale where log(N_H) + 12 = 12.
#
# Representative metals sufficient to evaluate f_α are included below.
# Omitted trace elements contribute negligibly to Z and do not affect f_α at the
# precision retained here.

using Test: @test
using StellarTracks.BaSTIv1: BaSTIv1Chemistry, alpha_mass_fraction

# ────────────────────────────────────────────────────────────────────────────
# Grevesse & Noels (1993) log-epsilon abundances: ε_i = log(N_i/N_H) + 12
# Columns: element => (atomic mass number A, log_epsilon)
# Note: the iron abundance log ε(Fe) = 7.67 is notably higher than in GS98 (7.50),
# which raises the non-alpha metal mass and thereby lowers f_α relative to GS98.
# Ar uses atomic mass 36 (the dominant stable alpha-process nucleus, ³⁶Ar) and is
# listed for completeness; it is NOT counted as an alpha element in BaSTIv1.
# ────────────────────────────────────────────────────────────────────────────
const _gn93 = (
    H  = ( 1, 12.00),
    He = ( 4, 10.99),
    C  = (12,  8.55),
    N  = (14,  7.97),
    O  = (16,  8.87),   # alpha
    Ne = (20,  8.09),   # alpha
    Na = (23,  6.33),
    Mg = (24,  7.58),   # alpha
    Al = (27,  6.47),
    Si = (28,  7.55),   # alpha
    P  = (31,  5.45),
    S  = (32,  7.21),   # alpha
    Ar = (36,  6.52),   # NOT alpha for BaSTIv1 (not modified in alpha-enhanced models)
    Ca = (40,  6.36),   # alpha
    Ti = (48,  5.02),   # alpha
    Fe = (56,  7.67),   # GN93 value; higher than GS98 (7.50)
    Ni = (58,  6.25),
)

const _alpha_elements_gn93 = (:O, :Ne, :Mg, :Si, :S, :Ca, :Ti)  # Ar excluded

# Mass contribution of each element: A_i · 10^(ε_i − 12)
_mass(A, ε) = A * exp10(ε - 12.0)

_contributions_gn93 = NamedTuple{keys(_gn93)}(
    _mass(A, ε) for (A, ε) in values(_gn93)
)

Z_total_gn93 = sum(v for (el, v) in pairs(_contributions_gn93) if el ∉ (:H, :He))
Z_alpha_gn93 = sum(v for (el, v) in pairs(_contributions_gn93) if el ∈ _alpha_elements_gn93)

# f_α = fraction of total metal mass contributed by alpha elements
f_alpha_gn93 = Z_alpha_gn93 / Z_total_gn93

# The stored value alpha_mass_fraction(BaSTIv1Chemistry()) = 0.6635 is f_alpha rounded
# to four decimal places.  We confirm the calculation is consistent with the stored value.
@test round(f_alpha_gn93; digits=4) == alpha_mass_fraction(BaSTIv1Chemistry())

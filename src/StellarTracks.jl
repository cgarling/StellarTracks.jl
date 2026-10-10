module StellarTracks

using ArgCheck: @argcheck
using TypedTables: Table

# For BCs.jl
using BolometricCorrections: AbstractBCTable, AbstractBCGrid, filternames, Mbol, logL, radius, surface_gravity, _parse_Mbol, MISTv2BCGrid
using BolometricCorrections.YBC: PARSECChemistry
import BolometricCorrections: AbstractChemicalMixture, X, X_phot, Y, Y_phot, Z, Z_phot, Y_p, MH, FeH, alphaFe, alpha_mass_fraction, chemistry, gridname
using StaticArrays: SVector
import Tables

# Top-level API definitions
# Track
""" Abstract supertype representing individual stellar tracks.
Different concrete implementations are used for different
stellar libraries, but all subtypes of `AbstractTrack` are guaranteed to support
the generic API described in the documentation.

The main way to get the properties of the
modeled star when it has a particular age, assuming `track` is a valid concrete instance,
is to call the track with the logarithmic age `track(log10(age))` where `age` is in years.
This will return a `NamedTuple` containing the available properties (e.g., logL, logTe)
of the modeled star at that age. Not all libraries will have the same properties available. """
abstract type AbstractTrack end 
Base.Broadcast.broadcastable(t::AbstractTrack) = Ref(t)
(track::AbstractTrack)(logAge::AbstractArray{<:Number}) = Table(track(la) for la in logAge)
# function (track::AbstractTrack)(logAge::AbstractArray{<:Number})
#     # result = track.itp(logAge)
#     # return Table(NamedTuple{(:logTe, :Mbol, :logg)}(i) for i in result)
#     # return Table(NamedTuple{(:logTe, :Mbol, :logg)}(track.itp(la)) for la in logAge)
#     # Conversion to Table is slightly slow ~40ns 
#     return Table(track(la) for la in logAge)
# end
"""
    gridname(::Type{<:AbstractTrack})

Returns a human-readable `String` identifier of the stellar evolution grid that tracks of the provided type are derived from. This should be provided for all track types.

    gridname(track::AbstractBCTable) = gridname(typeof(track))

Generic method provided to return the `gridname` for instances of concrete subtypes of `AbstractTrack`.
"""
gridname(t::AbstractTrack) = gridname(typeof(t))
"""
    Base.extrema(t::AbstractTrack)
Returns the minimum and maximum logarithmic age (`log10(age [yr])`) of the stellar model.
"""
function Base.extrema(t::AbstractTrack) end
"""
    Base.keys(t::AbstractTrack)
Returns a tuple of symbols that will be used to set the keys of the `NamedTuple`s returned by `t(logAge::Number)`.
"""
function Base.keys(t::AbstractTrack) end
"""
    mass(t::AbstractTrack)
Returns the initial stellar mass of the modeled star in solar masses. """
function mass(t::AbstractTrack) end
"""
    chemistry(t::AbstractTrack)
Returns an instance of a subtype of
[`AbstractChemicalMixture`](@extref BolometricCorrections.AbstractChemicalMixture)
valid for the provided track `t` that can be used in other methods.
"""
function chemistry(t::AbstractTrack) end
"""
    X(t::AbstractTrack)
Returns the hydrogen mass fraction of the modeled star. """
X(t::AbstractTrack) = X(chemistry(t), Z(t))
"""
    Y(t::AbstractTrack)
Returns the helium mass fraction of the modeled star. """
Y(t::AbstractTrack) = Y(chemistry(t), Z(t))
"""
    Z(t::AbstractTrack)
Returns the metal mass fraction of the modeled star. """
Z(t::AbstractTrack) = Z(chemistry(t), MH(t))
"""
    MH(t::AbstractTrack)
Returns the logarithmic metal abundance of the modeled star, defined as [M/H] = log(Z/X) - log(Z⊙/X⊙).
Generic method uses the [Salaris1993](@citet) equation to compute \\[M/H\\]
from [`FeH`](@ref) and [`alphaFe`](@ref). """
MH(t::AbstractTrack) = MH(chemistry(t), FeH(t), alphaFe(t)) # Calls generic from BolometricCorrections.jl
"""
    FeH(t::AbstractTrack)
Returns the iron abundance [Fe/H] of the modeled star. For scaled-solar models, this equals [`MH`](@ref).
Generic method uses the [Salaris1993](@citet) equation to compute \\[Fe/H\\]
from [`MH`](@ref) and [`alphaFe`](@ref)."""
FeH(t::AbstractTrack) = FeH(chemistry(t), MH(t), alphaFe(t)) # Calls generic from BolometricCorrections.jl
"""
    alphaFe(t::AbstractTrack)
Returns the α-element enhancement [α/Fe] of the modeled star. For scaled-solar models, this is zero. """
function alphaFe(t::AbstractTrack) end
"""
    post_rgb(t::AbstractTrack)
Returns `true` if the track includes post-RGB evolution, `false` otherwise. """
function post_rgb(t::AbstractTrack) end
###############################################
# TrackSet
""" Abstract supertype representing a set of stellar tracks
with common properties and supports
common operations like interpolating between tracks and
calculating isochrones. Specifically, it is assumed all individual
tracks in a track set have identical initial metallicity and chemical
composition. Subtypes of `AbstractTrackSet`
are guaranteed to support the generic API described in the documentation.

Concrete instances are callable with an initial stellar mass (in solar masses),
returning an interpolated track at the requested mass. """
abstract type AbstractTrackSet end
Base.Broadcast.broadcastable(ts::AbstractTrackSet) = Ref(ts)
# Canonical AMRs retain the high-mass side of physical folds, following Dotter (2016).
function _top_down_indices(masses::AbstractVector, ages::AbstractVector)
    order = sortperm(masses; rev=true)
    kept = Int[]
    for i in order
        if isempty(kept) || (masses[i] < masses[last(kept)] && ages[i] > ages[last(kept)])
            push!(kept, i)
        end
    end
    return kept
end
# Age-mass relation (AMR) at one EEP: log10(initial mass) linear in log age, on knots with masses decreasing as ages
# increase (e.g., from `_top_down_indices`). The interpolation weights depend only on the knot masses, so wherever two
# adjacent EEPs share knots, the isochrone initial mass cannot decrease from one to the next; PCHIP slopes depend on the
# ages and its curves can cross between shared knots. Log mass follows the near power-law lifetime-mass relation, so
# the accuracy is close to PCHIP's.
struct AMR{T <: AbstractVector, U <: AbstractVector}
    t::T # log age, increasing
    u::U # log10 initial mass, decreasing
end
_amr(m_ini, logAge) = AMR(logAge, log10.(m_ini))
# PCHIP interpolation of several properties sharing the knots `t` (initial masses, increasing), with SVector values `u`
# and slopes `du`
struct PCHIP{T <: AbstractVector, U <: AbstractVector}
    t::T
    u::U
    du::U
end
# Index j of the knot interval t[j] <= x < t[j+1] of increasing `t`, clamped to 1:length(t)-1; with `j`, found by
# walking from `j`, which is fast when it is close
_interval(t, x) = clamp(searchsortedlast(t, x), firstindex(t), lastindex(t) - 1)
function _interval(t, x, j)
    j = clamp(j, firstindex(t), lastindex(t) - 1)
    while j > firstindex(t) && t[j] > x; j -= 1; end
    while j < lastindex(t) - 1 && t[j+1] <= x; j += 1; end
    return j
end
_checkrange(t, x) = first(t) <= x <= last(t) || throw(DomainError(x, "outside the interpolation range ($(first(t)), $(last(t)))"))
# The arithmetic of the AMR and PCHIP evaluations is that of DataInterpolations' LinearInterpolation and
# CubicHermiteSpline, which these replaced
(a::AMR)(logAge) = (_checkrange(a.t, logAge); a(logAge, _interval(a.t, logAge)))
(a::AMR)(logAge, j) = a.u[j] + (a.u[j+1] - a.u[j]) / (a.t[j+1] - a.t[j]) * (logAge - a.t[j])
# Log age at log10 initial mass `lm` in the knot range; the inverse of the AMR, as linear interpolation on reversed knots
function _inverse(a::AMR, lm)
    t, u = a.t, a.u
    p = clamp(searchsortedfirst(u, lm; rev=true), firstindex(u) + 1, lastindex(u))
    return t[p] + (t[p-1] - t[p]) / (u[p-1] - u[p]) * (lm - u[p])
end
(p::PCHIP)(m) = (_checkrange(p.t, m); p(m, _interval(p.t, m)))
function (p::PCHIP)(m, k)
    t, u, du = p.t, p.u, p.du
    Δt = t[k+1] - t[k]
    c₁ = (u[k+1] - u[k] - du[k] * Δt) / Δt^2
    c₂ = (du[k+1] - du[k] - 2c₁ * Δt) / Δt^2
    Δt₀ = m - t[k]
    return u[k] + Δt₀ * du[k] + Δt₀^2 * (c₁ + (m - t[k+1]) * c₂)
end
# PCHIP slopes at knots `t` for values `u`: weighted harmonic mean of the neighboring secant slopes, or zero at local
# extrema, with the three-point end formula of pchiptx.m (Moler 2004, Numerical Computing with MATLAB, Section 3.6). Same
# arithmetic as DataInterpolations' internal `du_PCHIP`, which this replaced.
function _pchip_slopes(u, t)
    n = length(t)
    h = diff(t)
    δ = diff(u) ./ h
    function edge(h₁, h₂, δ₁, δ₂)
        d = ((2h₁ + h₂) * δ₁ - h₁ * δ₂) / (h₁ + h₂)
        sign(d) != sign(δ₁) && return zero(d)
        return sign(δ₁) != sign(δ₂) && abs(d) > 3abs(δ₁) ? 3δ₁ : d
    end
    du = similar(δ, n)
    du[1], du[n] = edge(h[1], h[2], δ[1], δ[2]), edge(h[n-1], h[n-2], δ[n-1], δ[n-2])
    for k in 2:n-1
        w₁, w₂ = 2h[k] + h[k-1], h[k] + 2h[k-1]
        du[k] = sign(δ[k-1]) == sign(δ[k]) != 0 ? (w₁ + w₂) / (w₁ / δ[k-1] + w₂ / δ[k]) : zero(eltype(δ))
    end
    return du
end
# The AMRs and property interpolants of all EEPs of a track set, stored flat: the knots of EEP i are `off[i]+1:off[i+1]`
# of one vector per field. `isochrone` visits every EEP, and reading them in sequence is several times faster than
# through separate objects. Indexing returns an `AMR` or `PCHIP` of views.
const _Segment{T} = SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true}
struct AMRSet{T, U} <: AbstractVector{AMR{_Segment{T}, _Segment{U}}}
    off::Vector{Int}
    t::Vector{T}
    u::Vector{U}
end
struct PCHIPSet{T, U} <: AbstractVector{PCHIP{_Segment{T}, _Segment{U}}}
    off::Vector{Int}
    t::Vector{T}
    u::Vector{U}
    du::Vector{U}
end
Base.size(s::Union{AMRSet, PCHIPSet}) = (length(s.off) - 1,)
Base.@propagate_inbounds Base.getindex(s::AMRSet, i::Int) = (r = s.off[i]+1:s.off[i+1]; AMR(view(s.t, r), view(s.u, r)))
Base.@propagate_inbounds Base.getindex(s::PCHIPSet, i::Int) = (r = s.off[i]+1:s.off[i+1]; PCHIP(view(s.t, r), view(s.u, r), view(s.du, r)))
# Vectors of the `AMR`s and `PCHIP`s of each EEP to flat storage
function _flatten(amrs, interps)
    amrs, interps = identity.(amrs), identity.(interps)
    off = [0; cumsum([length(a.t) for a in amrs])]
    cat(f, xs) = reduce(vcat, [getfield(x, f) for x in xs])
    return AMRSet(off, cat(:t, amrs), cat(:u, amrs)), PCHIPSet(off, cat(:t, interps), cat(:u, interps), cat(:du, interps))
end
# Indices of a longest strictly increasing subsequence of `x` (patience sorting, O(n log n))
function _increasing_subsequence(x::AbstractVector)
    tails = Int[] # tails[k] is the index of the smallest last element of an increasing subsequence of length k
    prev = zeros(Int, length(x))
    for i in eachindex(x)
        k = searchsortedfirst(tails, i; by=j -> x[j])
        prev[i] = k > 1 ? tails[k-1] : 0
        k > length(tails) ? push!(tails, i) : (tails[k] = i)
    end
    out = Int[]
    i = isempty(tails) ? 0 : last(tails)
    while i != 0
        push!(out, i)
        i = prev[i]
    end
    return reverse!(out)
end
# Whether the selected tracks (knots) of EEP indices i0 < i1 differ only because tracks end (selected at i0 and never
# after i1) or start (selected at i1 and never before i0)
function _coverage_change(interps, i0, i1)
    k0, k1 = interps[i0].t, interps[i1].t
    k0 == k1 && return false
    selected(x, r) = any(i -> insorted(x, interps[i].t), r)
    return all(x -> insorted(x, k1) || !selected(x, i1+1:lastindex(interps)), k0) &&
           all(x -> insorted(x, k0) || !selected(x, firstindex(interps):i0-1), k1)
end
# `_coverage_change` for each pair of neighboring EEPs, stored in the track set `properties` because it is too slow to
# evaluate for every isochrone
_coverage_flags(interps) = BitVector([_coverage_change(interps, i, i + 1) for i in firstindex(interps):lastindex(interps)-1])
# Where tracks end or start between the AMRs `a0` and `a1` of two EEPs, the isochrone mass step between them compares
# interpolations over different tracks and can be far from the physical step. At fixed age, a star of mass `m0` at the
# first EEP reaches the second after dt = t_1(m0) - t; this returns the ratio of the second mass to `m0` interpolated
# from the tracks present at both EEPs only if dt is negative or differs by more than a factor `K` from the time the
# two common tracks bracketing `m0` spend between the EEPs, and `nothing` otherwise.
function _common_step(a0, a1, logAge, m0; K=3)
    # Indices of the elements of sorted `xs` nearest to `x` on either side for which `keep` holds, or `nothing`
    function _bracket(keep, xs, x; rev=false)
        lo = searchsortedlast(xs, x; rev)
        lo == lastindex(xs) && xs[lo] == x && (lo -= 1)
        hi = lo + 1
        while lo >= firstindex(xs) && !keep(lo); lo -= 1; end
        while hi <= lastindex(xs) && !keep(hi); hi += 1; end
        return lo < firstindex(xs) || hi > lastindex(xs) ? nothing : (lo, hi)
    end
    _lerp(x, xs, ys, (lo, hi)) = ys[lo] + (x - xs[lo]) / (xs[hi] - xs[lo]) * (ys[hi] - ys[lo])
    in1(k) = insorted(a0.u[k], a1.u; rev=true) # knot k of `a0` is a track also selected at `a1`
    in0(k) = insorted(a1.u[k], a0.u; rev=true)
    lm = log10(m0)
    b, bm = _bracket(Returns(true), a1.u, lm; rev=true), _bracket(in1, a0.u, lm; rev=true)
    b0, b1 = _bracket(in1, a0.t, logAge), _bracket(in0, a1.t, logAge)
    (isnothing(b) || isnothing(bm) || isnothing(b0) || isnothing(b1)) && return nothing
    dt = exp10(_lerp(lm, a1.u, a1.t, b)) - exp10(logAge)
    d = map(k -> exp10(a1.t[searchsortedfirst(a1.u, a0.u[k]; rev=true)]) - exp10(a0.t[k]), bm)
    (dt <= 0 || dt > K * maximum(d) || dt < minimum(d) / K) || return nothing
    return exp10(_lerp(logAge, a1.t, a1.u, b1) - _lerp(logAge, a0.t, a0.u, b0))
end
# Indices into `ts.eeps` of the points of the isochrone of `ts` at `logAge`, the initial masses from the AMRs at which to
# evaluate properties, the initial masses to report, and the knot intervals of the property interpolants (`PCHIP`) at
# the first masses. Where tracks end or start between neighboring EEPs and the
# mass step between them is non-physical (`_common_step`), the reported masses replace that step with the step from
# the common tracks if this changes it by at most 1%, and shift all later masses by the same amount. Where the AMRs of
# neighboring EEPs use different knots the masses can still decrease along the EEPs (a fold); keeping the longest
# increasing subsequence drops the shorter side of each fold, so a small fold at the start of a long phase does not
# remove the phase.
function _isochrone_masses(ts::AbstractTrackSet, logAge::Number)
    idx, seg = Int[], Int[]
    masses = Vector{eltype(ts)}(undef, 0)
    j = 1
    for (i, amr) in enumerate(ts.AMRs)
        if first(amr.t) <= logAge <= last(amr.t) # Valid age range for the EEP
            # Neighboring EEPs have nearby knot intervals at a given age, so each search starts from the last
            j = _interval(amr.t, logAge, j)
            # The AMR and this EEP's property interpolant share knots in reverse order, so the clamp only absorbs
            # rounding in exp10 and the interval of the mass is at or next to that of the AMR
            itp = ts.interps[i]
            m = clamp(exp10(amr(logAge, j)), first(itp.t), last(itp.t))
            push!(idx, i)
            push!(masses, m)
            push!(seg, _interval(itp.t, m, lastindex(amr.t) - j))
        end
    end
    m_ini = copy(masses)
    offset = zero(eltype(m_ini))
    for j in firstindex(masses)+1:lastindex(masses)
        m_ini[j] = masses[j] + offset
        i0, i1 = idx[j-1], idx[j]
        (i1 == i0 + 1 ? ts.properties.coverage[i0] : _coverage_change(ts.interps, i0, i1)) || continue
        ratio = _common_step(ts.AMRs[i0], ts.AMRs[i1], logAge, masses[j-1])
        (isnothing(ratio) || abs(ratio * masses[j-1] / masses[j] - 1) > 0.01) && continue
        m_ini[j] = m_ini[j-1] * ratio
        offset = m_ini[j] - masses[j]
    end
    # Most isochrones are already strictly increasing
    issorted(m_ini; lt=<=) && return idx, masses, m_ini, seg
    keep = _increasing_subsequence(m_ini)
    return idx[keep], masses[keep], m_ini[keep], seg[keep]
end
# Masses for which at least one EEP has a valid property interpolation. Top-down selection can
# drop a mass from every EEP of a track set, so `unique(data.m_ini)` may advertise unusable masses.
_supported_masses(interps, m_ini) = filter(m -> any(itp -> first(itp.t) <= m <= last(itp.t), interps), unique(m_ini))
_pchip_svector(t, cols...) = PCHIP(t, SVector.(cols...), SVector.(map(c -> _pchip_slopes(c, t), cols)...))
# Generic function to interpolate a trackset to a new initial stellar mass
# This returns a NamedTuple with columns `(:logAge, names...)`, where `names` labels the
# components of the SVector-valued `ts.interps`; concrete subtypes should define methods to
# construct the correct AbstractTrack type from this NamedTuple.
function _generic_trackset_interp(ts::AbstractTrackSet, M::Number, names::NTuple{N, Symbol}) where {N}
    # Validate that mass is in range
    # throw(DomainError(M, "Requested mass $M is outside the valid range $(extrema(mass(ts))) for the track set."))
    m_min, m_max = extrema(mass(ts))
    @argcheck m_min <= M <= m_max "Requested mass $M is outside the valid range $(extrema(mass(ts))) for the track set."
    # Find EEP points where the requested mass is valid
    good_idx = findall(itp -> first(itp.t) <= M <= last(itp.t), ts.interps)
    # Convert to eltype(ts) so e.g. Float32 track sets return Float32 properties when evaluated at a Float64 mass
    results = [convert(SVector{N, eltype(ts)}, ts.interps[j](M)) for j in good_idx]
    # The clamp absorbs rounding in log10 at the ends of the mass range
    ages = [(amr = ts.AMRs[j]; _inverse(amr, clamp(log10(M), last(amr.u), first(amr.u)))) for j in good_idx]
    sortidx = sortperm(ages)
    sorted = results[sortidx]
    return NamedTuple{(:logAge, names...)}((ages[sortidx], ntuple(k -> getindex.(sorted, k), N)...))
end
"""
    gridname(::Type{<:AbstractTrackSet})

Returns a human-readable `String` identifier of the stellar evolution grid that track sets of the provided type are derived from. This should be provided for all track set types.

    gridname(trackset::AbstractTrackSet) = gridname(typeof(trackset))

Generic method provided to return the `gridname` for instances of concrete subtypes of `AbstractTrackSet`.
"""
gridname(ts::AbstractTrackSet) = gridname(typeof(ts))
"""
    mass(ts::AbstractTrackSet)
Returns the initial stellar masses (in solar masses) of the individual tracks contained in the track set. """
function mass(ts::AbstractTrackSet) end
"""
    chemistry(ts::AbstractTrack)
Returns an instance of a subtype of
[`AbstractChemicalMixture`](@extref BolometricCorrections.AbstractChemicalMixture)
valid for the provided track set `ts` that can be used in other methods.
"""
function chemistry(ts::AbstractTrackSet) end
"""
    X(ts::AbstractTrackSet)
Returns the common hydrogen mass fraction of the tracks contained in the track set. """
X(ts::AbstractTrackSet) = X.(Ref(chemistry(ts)), Z(ts))
"""
    Y(ts::AbstractTrackSet)
Returns the common helium mass fraction of the tracks contained in the track set. """
Y(ts::AbstractTrackSet) = Y.(Ref(chemistry(ts)), Z(ts))
"""
    Z(ts::AbstractTrackSet)
Returns the common metal mass fraction of the tracks contained in the track set. """
Z(ts::AbstractTrackSet) = Z.(Ref(chemistry(ts)), MH(ts))
"""
    MH(ts::AbstractTrackSet)
Returns the common logarithmic metal abundance of the tracks contained in the track set,
defined as [M/H] = log(Z/X) - log(Z⊙/X⊙).
Generic method uses the [Salaris1993](@citet) equation to compute \\[M/H\\]
from [`FeH`](@ref) and [`alphaFe`](@ref). """
MH(ts::AbstractTrackSet) = MH(chemistry(ts), FeH(ts), alphaFe(ts)) # Calls generic from BolometricCorrections.jl
"""
    FeH(ts::AbstractTrackSet)
Returns the common iron abundance [Fe/H] of the tracks contained in the track set.
For scaled-solar models, this equals [`MH`](@ref). Generic method uses the
[Salaris1993](@citet) equation to compute \\[Fe/H\\] from [`MH`](@ref) and
[`alphaFe`](@ref)."""
FeH(ts::AbstractTrackSet) = FeH(chemistry(ts), MH(ts), alphaFe(ts)) # Calls generic from BolometricCorrections.jl
"""
    alphaFe(ts::AbstractTrackSet)
Returns the common α-element enhancement [α/Fe] of the tracks contained in the track set. For scaled-solar models, this is zero. """
function alphaFe(ts::AbstractTrackSet) end
"""
    post_rgb(ts::AbstractTrackSet)
Returns `true` if the tracks in the track set include post-RGB evolution, `false` otherwise. """
function post_rgb(ts::AbstractTrackSet) end
"""
    isochrone(ts::AbstractTrackSet, logAge::Number)
Interpolates properties of the stellar tracks in the track set at
the requested logarithmic age (`logAge = log10(age [yr])`).
Returns a `NamedTuple` containing the properties; different libraries
may have different properties available. If `result = isochrone(...)`,
EEP points are generally available as `result.eep` and the corresponding
initial stellar masses are `result.m_ini`.
Isochrones contain only living stars (no stellar remnants), so the largest
`result.m_ini` is the initial mass of the most evolved living star; the exception is
a [`MISTv2TrackSet`](@ref StellarTracks.MIST.MISTv2TrackSet) created with `wd = true`.

```jldoctest
julia> using StellarTracks.MIST

julia> ts = MISTv1TrackSet(-1, 0); # Load set of MIST tracks with [M/H] = -1, vvcrit=0

julia> isochrone(ts, 10.0) isa NamedTuple
true
```
"""
function isochrone(ts::AbstractTrackSet, logAge::Number) end

###############################################
# TrackLibrary
""" Abstract supertype for loading the *full* set of all
stellar tracks available from a given library (e.g., PARSEC).
For some libraries (e.g., MIST), there are different model sets
for other stellar input parameters (e.g., rotation `vvcrit` in
MIST) that are not interpolated over -- these should be provided
as arguments to the concrete type constructors.

These are typically assembled as a collection of individual
track sets, with one track set constructed for each unique
set of stellar chemical compositions (i.e., [M/H]).
For more details, see the documentation for
[`AbstractTrackSet`](@ref StellarTracks.AbstractTrackSet).

Concrete subtypes must implement an isochrone interpolation method.
A generic method for the call signature
[`isochrone(tracklib::AbstractTrackLibrary, logAge::Number, mh::Number)`](@ref isochrone(::AbstractTrackLibrary, ::Number, ::Number))
is provided. This generic method should work as long as the concrete
subtype follows the [chemistry API](@ref chemistry_api) and has a field `tracklib.ts`
containing the individual tracksets that make up the library. The tracksets must support
isochrone interpolation via a method `isochrone(ts::<your_type_here>, logAge::Number)`.

Concrete instances must support interpolation of stellar tracks
at user-provided metallicities and initial stellar masses. A generic method
is implemented via the signature
`(tracklib::AbstractTrackLibrary)(mh::Number, M::Number)`. This method
interpolates the `tracklib` at logarithmic metallicity [M/H] = `mh` and stellar
initial mass `M`, returning an [`InterpolatedTrack`](@ref StellarTracks.InterpolatedTrack)
that can be called to evaluate the track."""
abstract type AbstractTrackLibrary end

"""
    gridname(::Type{<:AbstractTrackLibrary})

Returns a human-readable `String` identifier of the stellar evolution grid type provided. This should be provided for all track set types.

    gridname(tracklib::AbstractTrackLibrary) = gridname(typeof(tracklib))

Generic method provided to return the `gridname` for instances of concrete subtypes of `AbstractTrackLibrary`.
"""
gridname(tl::AbstractTrackLibrary) = gridname(typeof(tl))
"""
    chemistry(tl::AbstractTrackLibrary)
Returns an instance of a subtype of
[`AbstractChemicalMixture`](@extref BolometricCorrections.AbstractChemicalMixture)
valid for the provided track library `tl` that can be used in other methods.
"""
function chemistry(tl::AbstractTrackLibrary) end
"""
    X(tl::AbstractTrackLibrary)
Returns the hydrogen mass fractions of the track sets 
contained in the track library. """
X(tl::AbstractTrackLibrary) = X.(Ref(chemistry(tl)), Z(tl))
"""
    Y(tl::AbstractTrackLibrary)
Returns the helium mass fractions of the track sets
contained in the track library. """
Y(tl::AbstractTrackLibrary) = Y.(Ref(chemistry(tl)), Z(tl))
"""
    MH(tl::AbstractTrackLibrary)
Returns the metal mass fractions of the track sets
contained in the track library. """
Z(tl::AbstractTrackLibrary) = Z.(Ref(chemistry(tl)), MH(tl))
"""
    MH(tl::AbstractTrackLibrary)
Returns the logarithmic metal abundances of the track sets contained in the track library, 
defined as [M/H] = log(Z/X) - log(Z⊙/X⊙). Generic method uses the [Salaris1993](@citet) 
equation to compute \\[M/H\\] from [`FeH`](@ref) and [`alphaFe`](@ref). """
MH(tl::AbstractTrackLibrary) = MH.(Ref(chemistry(tl)), FeH(tl), alphaFe(tl)) # Calls generic from BolometricCorrections.jl
"""
    FeH(tl::AbstractTrackLibrary)
Returns the iron abundances [Fe/H] of the track sets contained in the track library.
For scaled-solar models, this equals [`MH`](@ref). Generic method uses the
[Salaris1993](@citet)  equation to compute \\[Fe/H\\] from [`MH`](@ref) and
[`alphaFe`](@ref)."""
FeH(tl::AbstractTrackLibrary) = FeH.(Ref(chemistry(tl)), MH(tl), alphaFe(tl)) # Calls generic from BolometricCorrections.jl
"""
    alphaFe(tl::AbstractTrackLibrary)
Returns the α-element enhancement [α/Fe] of the track sets in the track library. For scaled-solar models, this is zero. """
function alphaFe(tl::AbstractTrackLibrary) end
"""
    post_rgb(tl::AbstractTrackLibrary)
Returns `true` if any of the track sets
in the track library include post-RGB evolution, 
`false` otherwise. """
function post_rgb(tl::AbstractTrackLibrary) end

"""
    InterpolatedTrack(track0, track1, track1_prefac, track2_prefac) <: AbstractTrack

Type allowing for linear interpolation between two stellar tracks of identical 
mass but different metallicity (`track0` and `track1`). Simple linear 
interpolation between the two tracks is applied. The track can be evaluated 
at a particular logarithmic age by calling `(track::InterpolatedTrack)(logAge::Number)`. 

The constructor for this type is considered internal and users should not construct this type directly.
Rather, tracks of this type are constructed from track libraries by calling them with 
the signature `(tracklib::AbstractTrackLibrary)(mh::Number, M::Number)` where 
`mh` is the logarithmic metallicity [M/H] and `M` is the initial stellar mass in solar masses.

```jldoctest
julia> using StellarTracks.PARSEC: PARSECLibrary

julia> p = PARSECLibrary();

julia> track = p(-2.05, 1.05)
InterpolatedTrack with M_ini=1.05, MH=-2.05, Z=0.000138567, Y=0.248747, X=0.751115.

julia> track(9.0)
(logTe = 3.8820487347062302, Mbol = 3.7411721770340987, logg = 4.521853108813156, C_O = 0.0)
```
"""
struct InterpolatedTrack{A <: AbstractTrack, B <: Number} <: AbstractTrack
    track0::A
    track1::A
    track0_prefac::B
    track1_prefac::B
end
function (track::InterpolatedTrack)(logAge::Number)
    # Just let it fail if logAge errors for one of the two tracks
    result0 = track.track0(logAge)
    result1 = track.track1(logAge)
    result = values(result0) .* track.track0_prefac .+ values(result1) .* track.track1_prefac
    return NamedTuple{keys(track.track0)}(result)
end
gridname(::Type{<:InterpolatedTrack}) = throw(ArgumentError("InterpolatedTrack can be used with any track library, so a generic `gridname` is not defined. Call `gridname` on an concrete instance."))
gridname(track::InterpolatedTrack{A}) where {A} = gridname(A)
function Base.extrema(track::InterpolatedTrack)
    ext0 = extrema(track.track0)
    ext1 = extrema(track.track1)
    return (max(ext0[1], ext1[1]), min(ext0[2], ext1[2]))
end
Base.keys(track::InterpolatedTrack) = keys(track.track0)
Base.eltype(track::InterpolatedTrack) = promote_type(eltype(track.track0), eltype(track.track1))
mass(track::InterpolatedTrack) = mass(track.track0) == mass(track.track1) ? mass(track.track0) : mass(track.track0) * track.track0_prefac + mass(track.track1) * track.track1_prefac
chemistry(track::InterpolatedTrack) = chemistry(track.track0)
MH(track::InterpolatedTrack) = MH(track.track0) * track.track0_prefac + MH(track.track1) * track.track1_prefac
alphaFe(track::InterpolatedTrack) = alphaFe(track.track0) * track.track0_prefac + alphaFe(track.track1) * track.track1_prefac
post_rgb(track::InterpolatedTrack) = all(post_rgb, (track.track0, track.track1)) ? true : false
function Base.show(io::IO, mime::MIME"text/plain", t::InterpolatedTrack{A}) where A
    print(io, "InterpolatedTrack with M_ini=$(mass(t)), MH=$(round(MH(t); sigdigits=6)), Z=$(round(Z(t); sigdigits=6)), Y=$(round(Y(t); sigdigits=6)), X=$(round(X(t); sigdigits=6)).")
end
# function Base.show(io::IO, mime::MIME"text/plain", t::InterpolatedTrack{A}) where A
#     print(io, "Interpolation between two `$(split(string(A), "{")[1])`s with M_ini=$(mass(t)), MH=$(MH(t)), Z=$(Z(t)), Y=$(Y(t)), X=$(X(t)).")
# end

# Make AbstractTrackLibrary callable to construct interpolated track
function (tracklib::AbstractTrackLibrary)(mh::Number, M::Number)
    ts = tracklib.ts # vector of tracksets that make up tracklib
    mhvec = MH(tracklib)
    min_mh, max_mh = extrema(mhvec)
    if mh < min_mh || mh > max_mh
        throw(DomainError(mh, "Requested metallicity [M/H]=$mh is outside the valid range for the stellar track library of $(extrema(mhvec))."))
    end

    if !issorted(mhvec)
        idxs = sortperm(mhvec)
        ts = view(ts, idxs)
        mhvec = view(mhvec, idxs)
    end
    
    idx = searchsortedfirst(mhvec, mh)
    if mhvec[idx] ≈ mh # Requested mh matched in grid
        return ts[idx](M) # Assume a method (trackset)(M) exists
    else
        # Need to interpolate between two tracks with same mass, different metallicity
        track0 = ts[idx-1](M)
        track1 = ts[idx](M)
        mh0 = mhvec[idx-1]
        mh1 = mhvec[idx]
        return InterpolatedTrack(track0, track1, ((mh1-mh)/(mh1-mh0)), ((mh-mh0)/(mh1-mh0)))
    end
end

# This generic method was originally developed for MIST;
# should work for PARSEC as well after refactoring
"""
    isochrone(tracklib::AbstractTrackLibrary, logAge::Number, mh::Number)
Interpolates properties of the stellar tracks in the track library
at the requested logarithmic age (`logAge = log10(age [yr])`) and
metallicity [M/H] = `mh`. Returns a `NamedTuple` containing the
properties; different libraries may have different properties available.
If `result = isochrone(...)`, EEP points are generally available as
`result.eep` and the corresponding initial stellar masses are `result.m_ini`.
Between metallicity grid points, only EEPs present in both bracketing isochrones are returned.

```jldoctest
julia> using StellarTracks.MIST

julia> p = MISTv1Library(0.0); # Load the library of MIST tracks with vvcrit=0

julia> isochrone(p, 10.0, -1.65) isa NamedTuple
true
```
"""
function isochrone(tracklib::AbstractTrackLibrary, logAge::Number, mh::Number)
    ts = tracklib.ts
    mhvec = MH(tracklib)
    min_mh, max_mh = extrema(mhvec)
    if mh < min_mh || mh > max_mh
        throw(DomainError(mh, "Requested metallicity [M/H]=$mh is outside the valid range for the stellar track library of $(extrema(mhvec))."))
    end

    if !issorted(mhvec)
        idxs = sortperm(mhvec)
        ts = view(ts, idxs)
        mhvec = view(mhvec, idxs)
    end

    # searchsortedfirst returns the index of the first value in mhvec greater than or
    # equivalent to mh. If mh is greater than all values in mhvec, returns lastindex(mhvec) + 1.
    # We have already checked bounds so we know min_mh < mh < max_mh
    idx = searchsortedfirst(mhvec, mh)
    if mhvec[idx] ≈ mh # Requested mh matched in grid
        return isochrone(ts[idx], logAge)
    end

    # mh is valid, so need to interpolate isochrone as a function of [M/H]
    # According to Marigo2017, the interpolations (at least for BCs) in Z or [M/H]
    # are linear, so the BCs overall are computed by a
    # 3D linear interpolation in logTe x logg x [M/H] space.
    # VandenBerg2014 suggests linear interpolation is sufficient for Z or [M/H] interpolation,
    # although cubic interpolation was used in VandenBerg2012.
    # Most pre-computed grid seems more uniform in [M/H] than they are in Z, so I think it might
    # be a good idea to do linear interpolation in [M/H].
    
    # Evaluate isochrones on either side of intermediate point
    y0 = isochrone(ts[idx-1], logAge)
    y1 = isochrone(ts[idx], logAge)
    if length(first(y0)) == 0 || length(first(y1)) == 0
        throw(DomainError(logAge, "No valid EEPs were found for the requested logarithmic age $logAge."))
    end
    # Get indices into y0 and y1 that correspond to the overlapping EEP points
    y0_idxs = Vector{Int}(undef, 0)
    y1_idxs = similar(y0_idxs)
    good_eeps = similar(y0_idxs)
    for (y0_idx, eep) in enumerate(y0.eep)
        y1_idx = searchsortedfirst(y1.eep, eep)
        if y1_idx <= lastindex(y1.eep) && y1.eep[y1_idx] == eep
            push!(y0_idxs, y0_idx)
            push!(y1_idxs, y1_idx)
            push!(good_eeps, eep)
        end
    end
    # Get isochrone keys, removing EEP since that is fixed
    goodkeys = filter(Base.Fix1(!==, :eep), keys(y0))
    # Perform linear interpolation in _interp_kernel to establish a function
    # barrier to improve performance since some of the types of the variables
    # aren't known at runtime
    result = NamedTuple{goodkeys}(_interp_kernel(goodkeys, y0, y1, idx, y0_idxs, y1_idxs, mh, mhvec))
    # Concatenate interpolated result with valid EEP points
    return (eep = good_eeps, result...)
end
# This does linear interpolation for isochrone; function barrier for type stability
_interp_kernel(goodkeys, y0, y1, idx, y0_idxs, y1_idxs, x, xvec) =
    ((y0[key][y0_idxs] .* (xvec[idx] - x) .+ y1[key][y1_idxs] .* (x - xvec[idx-1])) ./ (xvec[idx] - xvec[idx-1]) for key in goodkeys)

###############################################
# Utilities

"""
    uniqueidx(v) = unique(Base.Fix1(getindex, v), eachindex(v))
Returns the indices of the first occurrences of unique elements in the array `v`. """
uniqueidx(v) = unique(Base.Fix1(getindex, v), eachindex(v)) # utility

###############################################
# Include files containing submodules for different track libraries
include(joinpath("parsec", "parsec.jl"))
using .PARSEC

include(joinpath("mist", "mist.jl"))
using .MIST

include(joinpath("basti", "v1", "basti_v1.jl"))
using .BaSTIv1

include(joinpath("basti", "v2", "basti_v2.jl"))
using .BaSTIv2

# Include bolometric correction-related functionality
include("BCs.jl")

include("precompile.jl")

# Common API exports
export mass, chemistry, X, Y, Z, MH, FeH, alphaFe, alpha_mass_fraction, post_rgb, isochrone, gridname
# Submodule exports
export PARSECLibrary, MISTv1Library, MISTLibrary, MISTv2Library, BaSTIv1Library, BaSTIv2Library

end

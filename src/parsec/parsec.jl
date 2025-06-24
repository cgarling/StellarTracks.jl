"""StellarTracks.PARSEC provides access to the PARSEC V1.2S stellar tracks."""
module PARSEC

# imports from parent module
using ..StellarTracks: AbstractChemicalMixture, AbstractTrack, AbstractTrackSet, AbstractTrackLibrary,
                       uniqueidx, _generic_trackset_interp
import ..StellarTracks: mass, post_rgb, isochrone
import ..StellarTracks: X, X_phot, Y, Y_phot, Z, Z_phot, MH, chemistry

import CSV
using CodecZlib: GzipDecompressorStream
using DataDeps: register, DataDep, @datadep_str, unpack
using DelimitedFiles: readdlm
using Glob: glob # file pattern matching
import JLD2 # for saving files in binary format
using ProgressMeter: @showprogress
# import Tables # For Tables.matrix conversion, "1" compat
import Tar
using TypedTables: Table
# import p7zip_jll: p7zip

# Imports for core module code
using ArgCheck: @argcheck
using StaticArrays: SVector
# import Interpolations: interpolate, Gridded, Linear, deduplicate_knots!
using DataInterpolations: AbstractInterpolation, AkimaInterpolation, LinearInterpolation, CubicSpline, PCHIPInterpolation, CubicHermiteSpline
Base.extrema(a::AbstractInterpolation) = extrema(a.t) # utility

# Number of EEP points per evolutionary phase
# starting from beginning of pre-main sequence through
# to thermally-pulsating AGB.
const eep_lengths = (PMS_BEG = 200,  # beginning of PMS to beginning of MS
                     MS_BEG = 200,   # beginning of MS to logTe min on MS
                     MS_TMIN = 200,  # logTe min on MS to MSTO 
                     MS_TO = 500,    # MSTO to TRGB
                     RG_TIP = 30,    # TRGB to beginning core He fusion
                     HE_BEG = 500,   # beginning core He fusion to end core He fusion
                     END_CHEB = 100, # end of core he fusion to beginning of TP-AGB
                     TPAGB_BEG = 200)
# Indices where the different phases begin
const eep_idxs = NamedTuple{keys(eep_lengths)}((1, (cumsum(values(eep_lengths)[begin:end-1]) .+ 1)...))
# HB files can have length 600 (HE_beg to END_CHEB) or 800 (including TPAGB as well)

# Column names for the stellar evolutionary tracks;
# units are log[yr], solar masses, log[K], bolometric mag, solar surface gravity (cgs)
# logAge in HB files stars from the ZAHB. Mbol and logg are calculated from parsec output as
# Mbol = 4.77 - 2.5 * log L
# log g = -10.616 + log mass + 4.0 * log Teff - log L
const track_header = SVector{6, String}("logAge", "mass", "logTe", "Mbol", "logg", "C_O")
const track_header_symbols = Tuple(Symbol.(i) for i in track_header)
# Which columns to actually keep after reading track file; used in Track and track_table
const select_columns = (:logAge, :logTe, :Mbol, :logg, :C_O)
# Matrix columns to keep when reading track in track_matrix
const keepcols = SVector(1,2,3,4,5,6)
const track_type = Float64 # Float type to use to represent values
"""Valid helium mass fractions (Y) for PARSECv1.2S."""
const ygrid = track_type[0.249, 0.249, 0.249, 0.25, 0.252, 0.256, 0.259, 0.263, 0.267, 0.273, 0.279, 0.284, 0.302, 0.321, 0.356]
"""Valid metal mass fractions (Z) for PARSECv1.2S."""
const zgrid = track_type[0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.014, 0.017, 0.02, 0.03, 0.04, 0.06]

##########################################################################
# PARSEC chemistry
"""
    PARSECChemistry()
Returns a singleton struct representing the PARSEC chemical mixture model.
We presently only include scaled-solar models. The solar protostellar chemical
mixture for PARSEC was calibrated to reproduce solar photospheric observations
via a forward modeling approach (see section 3 of [Bressan2012](@citet)). The
full solar calibration assumed for PARSEC is given in Table 3 of [Bressan2012](@citet).
The distribution of heavy metals is taken from [Grevesse1998](@citet) and [Caffau2011](@citet) (see section 4 of [Bressan2012](@citet)).

```jldoctest
julia> using StellarTracks.PARSEC: PARSECChemistry, X, Y, Z, X_phot, Y_phot, Z_phot, MH;

julia> chem = PARSECChemistry();

julia> X(chem) + Y(chem) + Z(chem) ≈ 1 # solar protostellar values
true

julia> X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1 # solar photospheric values
true

julia> MH(chem, Z(chem) * 0.1) ≈ -0.9400696788068212
true

julia> Z(chem, -0.9400696788068212) ≈ Z(chem) * 0.1
true
```
"""
struct PARSECChemistry <: AbstractChemicalMixture end
# For PARSEC, a choice can be made as to whether the initial solar
# chemical composition is taken to be the observed reference value
# i.e., Z⊙, Y⊙ in Table 3 of Bressan2012, or the photospheric abundances
# of the best-fit solar calibration model, i.e., Zs, Ys in Table 3. 
# For consistency with PARSEC's conversion between Z and [M/H], we will
# assume the observed reference values.

X(mix::PARSECChemistry) = 1 - Y(mix) - Z(mix) # 0.70226
X_phot(mix::PARSECChemistry) = 1 - Y_phot(mix) - Z_phot(mix)  # 0.73616
Y(::PARSECChemistry) = 0.28 # Y_initial in Table 3 of Bressan2012
Y_phot(::PARSECChemistry) = 0.2485  # Y⊙ in Table 3 of Bressan2012
# Y_phot(::PARSECChemistry) = 0.24787 # Y_S in Table 3 of Bressan2012
Z(::PARSECChemistry) = 0.01774 # Z_initial in Table 3 of Bressan2012
Z_phot(::PARSECChemistry) = 0.01524 # 0.01774 # Z⊙ in Table 3 of Bressan2012
# Z_phot(::PARSECChemistry) = 0.01597 # Z_S in Table 3 of Bressan2012
Y_p(::PARSECChemistry) = 0.2485

Y(mix::PARSECChemistry, Zval) = Y_p(mix) + 178//100 * Zval # γ = 1.78
# X generic
MH(mix::PARSECChemistry, Zval) = log10(Zval / X(mix, Zval)) - log10(Z_phot(mix) / X_phot(mix))
# MH(mix::PARSECChemistry, Zval) = log10(Zval / X(mix, Zval) / Z(mix) * X(mix))
function Z(mix::PARSECChemistry, MHval)
    # [M/H] = log(Z/X) - log(Z/X)☉ with Z☉ = solz
    # Z/X = exp10( [M/H] + log(Z/X)☉ )
    # X = 1 - Y - Z
    # Y ≈ Y_p + γ * Z for parsec (see Y(mix::PARSECChemistry, Zval) above)
    # so X ≈ 1 - (Y_p + γ * Z) - Z = 1 - Y_p - (1 + γ) * Z
    # Substitute into line 2,
    # Z / (1 - Y_p - (1 + γ) * Z) = exp10( [M/H] + log(Z/X)☉ )
    # Z = (1 - Y_p - (1 + γ) * Z) * exp10( [M/H] + log(Z/X)☉ )
    # let A = exp10( [M/H] + log(Z/X)☉ )
    # Z = (1 - Y_p) * A - (1 + γ) * Z * A
    # Z + (1 + γ) * Z * A = (1 - Y_p) * A
    # Z (1 + (1 + γ) * A) = (1 - Y_p) * A
    # Z = (1 - Y_p) * A / (1 + (1 + γ) * A)
    zoverx = exp10(MHval + log10(Z_phot(mix) / X_phot(mix)))
    γ = 178//100
    return (1 - Y_p(mix)) * zoverx / (1 + (1 + γ) * zoverx)
end

##########################################################################

# Data download, organization, file parsing, etc.
include("init.jl")

##########################################################################

"""
    PARSECTrack(zval::Number, mass::Number)
`PARSECTrack` implements the [`AbstractTrack`](@ref StellarTracks.AbstractTrack)
interface for the PARSEC stellar evolution library.

Note that due to the organization of the PARSEC data files, this method requires
constructing a [`PARSECTrackSet`](@ref) and is therefore
not efficient if your aim is to construct multiple tracks of the same metallicity `zval`. In this
case, you should construct a [`PARSECTrackSet`](@ref) and call it with the masses you want, e.g.,
`ts = PARSECTrackSet(0.0001); ts.([0.12, 0.15])`. 
```jldoctest
julia> track = StellarTracks.PARSEC.PARSECTrack(0.0001, 0.15)
PARSECTrack with M_ini=0.15, MH=-2.191722058538173, Z=0.0001, Y=0.248678, X=0.7512220000000001.

julia> track(7.0) # interpolate track at log10(age [yr]) = 7
(logTe = 3.6015066653099757, Mbol = 8.518315848633081, logg = 4.464972304683626, C_O = 0.0)
```
"""
struct PARSECTrack{A,B,C} <: AbstractTrack
    data::Table{A}
    itp::B
    properties::C
end
function PARSECTrack(filename::AbstractString) # Constructor from filename, raw parsec track files
    # Check file exists
    @argcheck isfile(filename)
    props = file_properties(filename)
    # Load data file into table
    data = CSV.read(filename, Table;
                    skipto = 2, header = Vector(track_header), types = track_type,
                    select = SVector(select_columns))
    # Construct interpolator as a function of proper age
    # itp = interpolate((exp10.(data.logAge),), [SVector(values(i)[2:end]) for i in data], Gridded(Linear()))
    itp = CubicSpline([SVector(values(i)[2:end]) for i in data], exp10.(data.logAge);
                      cache_parameters=true)
    return PARSECTrack(data, itp, (M = props.M, Z = props.Z, HB = props.HB))
end
# Constructor taking a subtable from one of the .jld2 reprocessed files
function PARSECTrack(data::Table, zval::Number, mass::Number)
    return PARSECTrack(data, CubicSpline([SVector(values(i)[2:end]) for i in data], exp10.(data.logAge)),
                       (M = mass, Z = zval, HB = length(data) > eep_idxs.RG_TIP))
end
# Constructor taking Z value, initial stellar mass, loads Table, calls above method
function PARSECTrack(zval::Number, mass::Number)
    # For PARSEC, individual tracks are not saved, so we need to load a trackset
    data = PARSECTrackSet(zval)(mass).data
    return PARSECTrack(data, zval, mass) # Method above
end
# Make Track callable with logAge to get logTe, Mbol, and logg as a NamedTuple
Base.keys(::PARSECTrack) = select_columns[2:end]
function (track::PARSECTrack)(logAge::Number)
    result = track.itp(exp10(logAge))
    return NamedTuple{keys(track)}(result)
end
Base.extrema(t::PARSECTrack) = log10.(extrema(t.itp.t))
mass(t::PARSECTrack) = t.properties.M
chemistry(::PARSECTrack) = PARSECChemistry()
MH(t::PARSECTrack) = MH(chemistry(t), Z(t))
Z(t::PARSECTrack) = t.properties.Z
Y(t::PARSECTrack) = Y(chemistry(t), Z(t))
X(t::PARSECTrack) = 1 - Y(t) - Z(t)
post_rgb(t::PARSECTrack) = t.properties.HB #hashb()
Base.eltype(t::PARSECTrack) = typeof(t.properties.Z)
function Base.show(io::IO, mime::MIME"text/plain", t::PARSECTrack)
    print(io, "PARSECTrack with M_ini=$(mass(t)), MH=$(MH(t)), Z=$(Z(t)), Y=$(Y(t)), X=$(X(t)).")
end

##########################################################################

"""
    PARSECTrackSet(zval::Number)
`PARSECTrackSet` implements the [`AbstractTrackSet`](@ref StellarTracks.AbstractTrackSet)
interface for the PARSEC stellar evolution library.
```jldoctest
julia> ts = StellarTracks.PARSEC.PARSECTrackSet(0.0001)
TrackSet with Y=0.248678, Z=0.0001, 1930 EEPs and 104 initial stellar mass points.

julia> ts(1.01) # Interpolate track at new initial mass
PARSECTrack with M_ini=1.01, MH=-2.191722058538173, Z=0.0001, Y=0.248678, X=0.7512220000000001.

julia> isochrone(ts, 10.0) isa NamedTuple # Interpolate isochrone at `log10(age [yr]) = 10`
true
```
"""
struct PARSECTrackSet{A <: AbstractVector{<:Integer},
                      B <: AbstractVector{<:AbstractInterpolation}, # AbstractInterpolation{T}
                      # C <: AbstractVector{<:AbstractVector{<:AbstractInterpolation}},
                      C,
                      D} <: AbstractTrackSet
    eeps::A # EEP points; vector of integers
    AMRs::B # age-mass relations; vector of interpolators, same length as eeps
    interps::C
    properties::D
end
function PARSECTrackSet(data::Table, Z::Number)
    # Now that we have an orderly data matrix, we need to construct the age-mass relations
    # age (dependent) is a monotonic function of m_ini (independent)
    eeps = sort(unique(data.eep))
    # @time amrs = [begin
    #             tmpdata = data[data.eep .== eep]
    #             # `data` is already sorted by mass, now sort by age
    #             # which will be the independent variable in the interpolation
    #             idxs = sortperm(tmpdata.logAge)
    #             tmpdata = tmpdata[idxs]
    #             # Now keep only entries that define a unique mapping between age and m_ini
    #             tmpdata = tmpdata[uniqueidx(tmpdata.logAge)]
    #             # Bergbusch2001 and VanenBerg2012 use cubic / akima interpolation
    #             # because their age-mass relations are monotonic at fixed EEP.
    #             # These relations are not precisely monotonic in PARSEC so we are getting
    #             # some overshoot and non-physical behavior. try linear interp for now
    #             # CubicSpline([SVector(values(i)[2:5]) for i in tmpdata], exp10.(tmpdata.logAge))
    #             # CubicSpline(tmpdata.m_ini, exp10.(tmpdata.logAge))
    #             # CubicSpline(tmpdata.m_ini, tmpdata.logAge)
    #             # AkimaInterpolation(tmpdata.m_ini, tmpdata.logAge)
    #             # AkimaInterpolation(tmpdata.m_ini, exp10.(tmpdata.logAge))
    #             # LinearInterpolation(tmpdata.m_ini, tmpdata.logAge)
    #             # LinearInterpolation(tmpdata.m_ini, exp10.(tmpdata.logAge))
    #             # PCHIP nice because it never overshoots the data
    #             PCHIPInterpolation(tmpdata.m_ini, tmpdata.logAge)
    #             # [[SVector(values(i)[2:4]) for i in tmpdata], exp10.(tmpdata.logAge)]
    #         end
    #         for eep in eeps]
    # Now construct interpolator for L, T, etc. as a function of mass
    # pinterps = [begin
    #                 tmpdata = data[data.eep .== eep]
    #                 # Sort by initial mass
    #                 idxs = sortperm(tmpdata.m_ini)
    #                 tmpdata = tmpdata[idxs]
    #                 CubicSpline([SVector(values(i)[[2,3,4,5,6,7]]) for i in tmpdata],
    #                             tmpdata.m_ini)
    #             end
    #             for eep in eeps]
    # Maybe best to just do an interpolator for each property eh?
    # We should wrap this into a single loop
    itp_type1 = CubicHermiteSpline{Vector{track_type},
                                   Vector{track_type},
                                   Vector{track_type},
                                   Vector{track_type},
                                   Vector{track_type},
                                   track_type}
    itp_type2 = LinearInterpolation{Vector{track_type},
                                    Vector{track_type},
                                    Vector{track_type},
                                    Vector{track_type},
                                    track_type}
    itp_type3 = AkimaInterpolation{Vector{track_type}, Vector{track_type}, Vector{track_type}, Vector{Float64}, Vector{Float64}, Vector{Float64}, track_type}

    amrs = Vector{itp_type1}(undef, length(eeps))
    logte = Vector{itp_type1}(undef, length(eeps))
    mbol = similar(logte)
    logg = similar(logte)
    c_o = similar(logte)
    Threads.@threads for i in eachindex(eeps)
        eep = eeps[i]
        # tmpdata = data[data.eep .== eep]
        tmpdata = data[findall(Base.Fix1(==, eep), data.eep)] # Performance optimization
        # Sort by age, which will be the independent variable in the AMR interpolation
        idxs = sortperm(tmpdata.logAge)
        tmpdata = tmpdata[idxs]
        # Now keep only entries that define a unique mapping between age and m_ini
        tmpdata = tmpdata[uniqueidx(tmpdata.logAge)]
        
        # We'll start the interpolation at the most massive model and enforce a monotonic
        # decrease in stellar initial mass with increasing age at fixed EEP
        _, p1 = findmax(tmpdata.m_ini)
        idxs = p1:lastindex(tmpdata)
        goodidxs = diff(tmpdata.m_ini[idxs]) .< 0
        goodidxs = vcat(true, goodidxs) # add true for first element as well
        tmpdata = tmpdata[idxs[goodidxs]]

        # PCHIPInterpolation is a type of CubicHermiteSpline
        amrs[i] = PCHIPInterpolation(tmpdata.m_ini, tmpdata.logAge)
        # amrs[i] = PCHIPInterpolation(tmpdata.m_ini[idxs[goodidxs]], exp10.(tmpdata.logAge[idxs[goodidxs]]))
        # amrs[i] = AkimaInterpolation(tmpdata.m_ini[idxs[goodidxs]], exp10.(tmpdata.logAge[idxs[goodidxs]]))
        # amrs[i] = AkimaInterpolation(tmpdata.m_ini, exp10.(tmpdata.logAge))
        # amrs[i] = AkimaInterpolation(tmpdata.m_ini, tmpdata.logAge)
        # TBH using linear interpolation for the AMRs does not seem to be that bad even
        # amrs[i] = LinearInterpolation(tmpdata.m_ini[idxs[goodidxs]], exp10.(tmpdata.logAge[idxs[goodidxs]]))
        # amrs[i] = LinearInterpolation(tmpdata.m_ini, tmpdata.logAge)
        # amrs[i] = LinearInterpolation(tmpdata.m_ini, exp10.(tmpdata.logAge))
        # Sort by initial stellar mass for defining the other interpolations
        idxs = sortperm(tmpdata.m_ini)
        tmpdata = tmpdata[idxs]
        # Now interpolate mass against logte, mbol, logg, c_o
        logte[i] = PCHIPInterpolation(tmpdata.logTe, tmpdata.m_ini)
        mbol[i] = PCHIPInterpolation(tmpdata.Mbol, tmpdata.m_ini)
        logg[i] = PCHIPInterpolation(tmpdata.logg, tmpdata.m_ini)
        c_o[i] = PCHIPInterpolation(tmpdata.C_O, tmpdata.m_ini)
        # logte[i] = LinearInterpolation(tmpdata.logTe, tmpdata.m_ini)
        # mbol[i] = LinearInterpolation(tmpdata.Mbol, tmpdata.m_ini)
        # logg[i] = LinearInterpolation(tmpdata.logg, tmpdata.m_ini)
        # c_o[i] = LinearInterpolation(tmpdata.C_O, tmpdata.m_ini)
    end
    
    return PARSECTrackSet(eeps, amrs, (logTe = logte, Mbol = mbol, logg = logg, C_O = c_o),
                          (Z = Z, masses = unique(data.m_ini)))
end
function PARSECTrackSet(zval::Number)
    idx = findfirst(≈(zval), zgrid) # Validate against zgrid
    if isnothing(idx)
        throw(ArgumentError("Provided `zval` argument $zval to `PARSECTrackSet` is invalid; available metal mass fractions are $zvals. For metallicity interpolation, use `PARSECLibrary`."))
    end
    fname = "Z" * string(zgrid[idx]) * "Y" * string(ygrid[idx]) * ".jld2"
    dd = @datadep_str(joinpath("PARSECv1.2S", fname))
    table = JLD2.load_object(dd)
    return PARSECTrackSet(table, zval)
end
(ts::PARSECTrackSet)(M::Number) = PARSECTrack(Table(_generic_trackset_interp(ts, M)), Z(ts), M)
mass(ts::PARSECTrackSet) = ts.properties.masses
chemistry(::PARSECTrackSet) = PARSECChemistry()
Z(ts::PARSECTrackSet) = ts.properties.Z
Y(ts::PARSECTrackSet) = Y(chemistry(ts), Z(ts))
X(ts::PARSECTrackSet) = 1 - Y(ts) - Z(ts)
MH(ts::PARSECTrackSet) = MH(chemistry(ts), Z(ts))
post_rgb(ts::PARSECTrackSet) = ts.eeps[end] > eep_idxs.RG_TIP
Base.eltype(ts::PARSECTrackSet) = typeof(ts.properties.Z)
function Base.show(io::IO, mime::MIME"text/plain", ts::PARSECTrackSet)
    print(io, "TrackSet with Y=$(Y(ts)), Z=$(Z(ts)), $(length(ts.AMRs)) EEPs and $(length(mass(ts))) initial stellar mass points.")
end
function isochrone(ts::PARSECTrackSet, logAge::Number) # 800 μs
    eeps = Vector{Int}(undef, 0)
    track_extrema = extrema(mass(ts))
    # interp_masses = Vector{eltype(first(ts.AMRs).u)}(undef, 0)
    interp_masses = Vector{eltype(ts)}(undef, 0)
    logte = similar(interp_masses)
    mbol = similar(interp_masses)
    logg = similar(interp_masses)
    c_o = similar(interp_masses)
    # data = Vector{eltype(first(ts.AMRs))}(undef, 0)
    for (i, amr) in enumerate(ts.AMRs)
        drange = extrema(amr) # Get valid age range for the EEP
        if logAge >= first(drange) && logAge <= last(drange)
        # if age >= first(drange) && age <= last(drange)
            imass = amr(logAge) # 100 ns
            # imass = amr(age)
            # Due to non-linear interpolation of mass(age), it is possible for the
            # interpolated masses to be outside the range valid for the tracks
            # Check here and do not write if outside range. Also require monotonically
            # increasing initial masses (gets rid of some non-monotonic behavior in the
            # underlying EEP tracks). 
            if imass >= first(track_extrema) && imass <= last(track_extrema) # && (length(interp_masses) == 0 || imass > last(interp_masses))
                mboli = ts.interps.Mbol[i](imass) # 120 ns
                # Enforce monotonically increasing Mbol along the MS (which ends at eep_idxs[4])
                # We can't be sure if last point was overly bright or current point is overly faint,
                # so we delete last point and continue so this point isn't output
                # Maybe not necessary, revisit later
                if length(mbol) > 0 && i < eep_idxs[4] && mboli > last(mbol)
                    li = lastindex(mbol)
                    deleteat!(eeps, li)
                    deleteat!(interp_masses, li)
                    deleteat!(logte, li)
                    deleteat!(mbol, li)
                    deleteat!(logg, li)
                    deleteat!(c_o, li)
                    continue
                end
                push!(eeps, i)
                push!(interp_masses, imass)
                push!(mbol, mboli)
                push!(logte, ts.interps.logTe[i](imass))
                push!(logg, ts.interps.logg[i](imass))
                push!(c_o, ts.interps.C_O[i](imass))
            end
        end
    end
    return (eep = eeps, m_ini = interp_masses, logTe = logte, Mbol = mbol, logg = logg, C_O = c_o)
end

##########################################################################

"""
    PARSECLibrary()
`PARSECLibrary` implements the [`AbstractTrackLibrary`](@ref StellarTracks.AbstractTrackLibrary)
interface for the PARSEC stellar evolution library. If you construct an instance as
`p = PARSECLibrary()`, it is callable as `p(mh::Number, M::Number)` which returns 
an [`InterpolatedTrack`](@ref StellarTracks.InterpolatedTrack)
that interpolates between tracks to a specific metallicity ([M/H]) and initial stellar mass (`M`).

This type also supports isochrone construction
(see [isochrone](@ref StellarTracks.isochrone(::StellarTracks.PARSEC.PARSECLibrary, ::Number, ::Number))).

# Examples
```jldoctest
julia> p = PARSECLibrary()
Structure of interpolants for PARSEC v1.2S library of stellar tracks. Valid range of metal mass fraction Z is (0.0001, 0.06).

julia> isochrone(p, 10.05, -0.76) isa NamedTuple
true

julia> p(-2.05, 1.05)
InterpolatedTrack with M_ini=1.05, MH=-2.05, Z=0.00013856708164357998, Y=0.24874664940532557, X=0.7511147835130308.
```
"""
struct PARSECLibrary{A} <: AbstractTrackLibrary
    ts::A # Vector of `TrackSet`s
end
chemistry(::PARSECLibrary) = PARSECChemistry()
Z(p::PARSECLibrary) = zgrid
Y(p::PARSECLibrary) = Y.(chemistry(p), Z(p))
X(p::PARSECLibrary) = 1 .- Y(p) .- Z(p)
# MH(p::PARSECLibrary) = PARSEC_MH.(Z(p)) # MH.(Z(p), Y(p))
MH(tl::PARSECLibrary) = MH.(chemistry(tl), Z(tl))
post_rgb(t::PARSECLibrary) = true
Base.eltype(p::PARSECLibrary) = typeof(first(Z(p)))
Base.Broadcast.broadcastable(p::PARSECLibrary) = Ref(p)
function Base.show(io::IO, mime::MIME"text/plain", p::PARSECLibrary)
    print(io, "Structure of interpolants for PARSEC v1.2S library of stellar tracks. Valid range of metal mass fraction Z is $(extrema(Z(p))).")
end
function PARSECLibrary()
    # Processing into TrackSet structures takes additional time
    # set_files = glob("Z*.jld2", base_dir)
    # filestems = [splitext(basename(fi))[1] for fi in set_files]
    # Z = [parse(Float64, split(split(fi, 'Z')[2], 'Y')[1]) for fi in filestems]
    # Y = [parse(Float64, split(fi, 'Y')[2]) for fi in filestems]
    # # Sort according to Z; isochrone(p::PARSEC...) depends on this
    # idxs = sortperm(Z)
    # Z .= Z[idxs]
    # Y .= Y[idxs]
    # set_files .= set_files[idxs]
    # # Make vector of tracksets
    # ts = [PARSECTrackSet(JLD2.load_object(fname), Z[i], Y[i]) for (i, fname) in enumerate(set_files)]
    # return PARSECLibrary(ts, Z, Y)
    return PARSECLibrary(PARSECTrackSet.(zgrid))
end
"""
    isochrone(p::PARSECLibrary, logAge::Number, mh::Number)
Interpolates properties of the stellar tracks in the library at the requested logarithmic age (`logAge = log10(age [yr])`) and logarithmic metallicity `mh`. Returns a `NamedTuple` containing the properties listed below:
 - `eep`: Equivalent evolutionary points
 - `m_ini`: Initial stellar masses, in units of solar masses.
 - `logTe`: Base-10 logarithm of the effective temperature [K] of the stellar model.
 - `Mbol`: Bolometric luminosity of the stellar model.
 - `logg`: Surface gravity of the stellar model calculated as `-10.616 + log10(mass) + 4 * logTe - (4.77 - Mbol) / 2.5`.
 - `C_O`: Photospheric C/O ratio (the ZAMS value is used before the TP-AGB).
"""
isochrone(p::PARSECLibrary, logAge::Number, mh::Number) # Falls back to generic method in StellarTracks.jl

#################################################################################

export PARSECTrack, PARSECTrackSet, PARSECLibrary, PARSECChemistry # Unique module exports
export mass, chemistry, X, Y, Z, MH, post_rgb, isochrone # Export generic API methods

#################################################################################

# padova_tracks does something to interpolate between the base track and the HB tracks
# by doing a linear interpolation between the two with a number of points equal to the
# number used for the HE_BEG stage which is 500. This feels like overkill to me ...
# in the parsec isochrones it's just a step function from the TRGB to the ZAHB
# as the helium flash is extremely short ~10,000 years.
# https://www.atnf.csiro.au/outreach/education/senior/astrophysics/stellarevolution_postmain.html#postmainheflash
# Technically the ZAHB is formed by the spread in mass-loss rates of stars transitioning
# off of the RGB, so maybe mass is the main factor there, not age?

# PARSEC webform will let you request Zini < 0.0001, but the output will always have
# MH ≈ -2.2 even when Zini is much lower. The theoretical quantities from the tracks
# e.g., logTe and logL are the same for Zini < 0.0001 and Zini=0.0001, indicating
# that the underlying stellar tracks do not change for Zini < 0.0001. However,
# the isochrone magnitudes *do* change, indicating that the webform may be applying
# the lower Zini to the bolometric corrections while using Zini=0.0001 for the stellar track.
# This would make sense if you expected little change in the stellar track as a function of
# Z for Zini < 0.0001.

end # module

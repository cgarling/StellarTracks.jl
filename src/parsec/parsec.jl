"""StellarTracks.PARSEC provides access to the PARSEC V1.2S stellar tracks."""
module PARSEC

# imports from parent module
using ..StellarTracks: AbstractChemicalMixture, AbstractTrack, AbstractTrackSet, AbstractTrackLibrary, uniqueidx
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
const select_columns = SVector(:logAge, :mass, :logTe, :Mbol, :logg, :C_O)
# Matrix columns to keep when reading track in track_matrix
const keepcols = SVector(1,2,3,4,5,6)
const track_type = Float64 # Float type to use to represent values

##########################################################################
# PARSEC chemistry
"""
    PARSECChemistry()
Returns a singleton struct representing the PARSEC chemical mixture model.
We presently only include scaled-solar models. The solar protostellar chemical
mixture for PARSEC was calibrated to reproduce solar photospheric observations
via a forward modeling approach (see section 3 of [Bressan2012](@citet)). The
full solar calibration assumed for PARSEC is given in Table 3 of [Bressan2012](@citet).
The distribution of heavy metals is taken from [Grevesse1998](@citet) and [Caffau2011](@citet)
(see section 4 of [Bressan2012](@citet)).

```jldoctest
julia> using StellarTracks.PARSEC: PARSECChemistry, X, Y, Z, X_phot, Y_phot, Z_phot, MH;

julia> chem = PARSECChemistry();

julia> X(chem) + Y(chem) + Z(chem) ≈ 1 # solar protostellar values
true

julia> X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1 # solar photospheric values
true

julia> MH(chem, Z(chem) * 0.1) ≈ -1.0265716016323736
true

julia> Z(chem, -1.0265716016323736) ≈ Z(chem) * 0.1
true
```
"""
struct PARSECChemistry <: AbstractChemicalMixture end
X(::PARSECChemistry) = 0.70226
X_phot(::PARSECChemistry) = 0.73616
Y(::PARSECChemistry) = 0.28 # Y_initial in Table 3 of Bressan2012
Y_phot(::PARSECChemistry) = 0.24787 # Y_S in Table 3 of Bressan2012
Y_p(::PARSECChemistry) = 0.2485
Z(::PARSECChemistry) = 0.01774 # Z_initial in Table 3 of Bressan2012
Z_phot(::PARSECChemistry) = 0.01597 # Z_S in Table 3 of Bressan2012

Y(mix::PARSECChemistry, Zval) = Y_p(mix) + 1.78 * Zval # γ = 1.78
# X generic
MH(mix::PARSECChemistry, Zval) = log10(Zval / X(mix, Zval)) - log10(Z(mix) / X(mix))
# MH(mix::PARSECChemistry, Zval) = log10(Zval / X(mix, Zval) / Z(mix) * X(mix))
"""
    MH_canon(mix::PARSECChemistry, Zval)
Returns the [M/H] value for the provided protostellar value of `Z`, where the solar abundances
are taken to be *photospheric* rather than protostellar, while the provided `Zval` is assumed
to be protostellar. This matches the convention of Table 4
of [Bressan2012](@citet) and the CMD webform for PARSEC.

```jldoctest
julia> using StellarTracks.PARSEC: PARSECChemistry, MH_canon;

julia> MH_canon(PARSECChemistry(), 0.0005) ≈ -1.4921252963659897
true
```
"""
MH_canon(mix::PARSECChemistry, Zval) = log10(Zval / X(mix, Zval)) - log10(0.0207)
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
    zoverx = exp10(MHval + log10(Z(mix) / X(mix)))
    γ = 1.78
    return (1 - Y_p(mix)) * zoverx / (1 + (1 + γ) * zoverx)
end
"""
    Z_canon(mix::PARSECChemistry, MHval)
Returns the *protostellar* metal mass fraction corresponding to the provided value
of [M/H] where the solar abundances are taken to be *photospheric* rather
than protostellar. This matches the convention of Table 4
of [Bressan2012](@citet) and the CMD webform for PARSEC.

```jldoctest
julia> using StellarTracks.PARSEC: PARSECChemistry, Z_canon;

julia> Z_canon(PARSECChemistry(), -1.4921252963659897) ≈ 0.0005
true
```
"""
function Z_canon(mix::PARSECChemistry, MHval)
    zoverx = exp10(MHval + log10(0.0207))
    γ = 1.78
    return (1 - Y_p(mix)) * zoverx / (1 + (1 + γ) * zoverx)
end

##########################################################################

# Data download, organization, file parsing, etc.
include("init.jl")

##########################################################################

""" `PARSECTrack` implements the [`AbstractTrack`](@ref StellarTracks.AbstractTrack)
interface for the PARSEC stellar evolution library. """
struct PARSECTrack{A,B,C} <: AbstractTrack
    filename::String
    data::Table{A}
    itp::B
    properties::C
end
function PARSECTrack(filename::AbstractString) # Constructor from filename
    # Check file exists
    @argcheck isfile(filename)
    props = file_properties(filename)
    # Load data file into table
    data = CSV.read(filename, Table;
                    skipto = 2, header = Vector(track_header), types = track_type,
                    select = select_columns)
    # Construct interpolator as a function of proper age
    # itp = interpolate((exp10.(data.logAge),), [SVector(values(i)[2:end]) for i in data], Gridded(Linear()))
    itp = CubicSpline([SVector(values(i)[2:end]) for i in data], exp10.(data.logAge))
    return PARSECTrack(convert(String, filename), data, itp, props)
end
# Make Track callable with logAge to get logTe, Mbol, and logg as a NamedTuple
function (track::PARSECTrack)(logAge::Number)
    # result = track.itp(logAge)
    result = track.itp(exp10(logAge))
    # return NamedTuple{Tuple(Symbol(i) for i in track_header[3:5])}(result)
    # return NamedTuple{(:logTe, :Mbol, :logg)}(result)
    return NamedTuple{Tuple(select_columns)[2:end]}(result)
end
function (track::PARSECTrack)(logAge::AbstractArray{<:Number})
    # result = track.itp(logAge)
    # return Table(NamedTuple{(:logTe, :Mbol, :logg)}(i) for i in result)
    # return Table(NamedTuple{(:logTe, :Mbol, :logg)}(track.itp(la)) for la in logAge)
    # Conversion to Table is slightly slow ~40ns 
    return Table(track(la) for la in logAge)
end
Base.extrema(t::PARSECTrack) = log10.(extrema(t.itp.t))
mass(t::PARSECTrack) = t.properties.M
Z(t::PARSECTrack) = t.properties.Z
Y(t::PARSECTrack) = t.properties.Y
X(t::PARSECTrack) = 1 - Y(t) - Z(t)
chemistry(::PARSECTrack) = PARSECChemistry()
MH(t::PARSECTrack) = MH(chemistry(t), Z(t))
post_rgb(t::PARSECTrack) = t.properties.HB #hashb()
Base.eltype(t::PARSECTrack) = typeof(t.properties.Z)

##########################################################################

""" `PARSECTrackSet` implements the [`AbstractTrackSet`](@ref StellarTracks.AbstractTrackSet)
interface for the PARSEC stellar evolution library. """
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
function PARSECTrackSet(data::Table, Z::Number, Y::Number)
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
    
    # logte = PCHIP
    # return pinterps
    # return TrackSet(eeps, amrs, [logte, mbol, logg, c_o], (Z = Z, Y = Y, masses = ms_props.M))
    return PARSECTrackSet(eeps, amrs, (logTe = logte, Mbol = mbol, logg = logg, C_O = c_o),
                          (Z = Z, Y = Y, masses = unique(data.m_ini)))
end
mass(ts::PARSECTrackSet) = ts.properties.masses
Z(ts::PARSECTrackSet) = ts.properties.Z
Y(ts::PARSECTrackSet) = ts.properties.Y
X(ts::PARSECTrackSet) = 1 - Y(ts) - Z(ts)
chemistry(::PARSECTrackSet) = PARSECChemistry()
MH(ts::PARSECTrackSet) = MH(chemistry(ts), Z(ts))
Base.eltype(ts::PARSECTrackSet) = typeof(ts.properties.Z)
function Base.show(io::IO, mime::MIME"text/plain", ts::PARSECTrackSet)
    print(io, "TrackSet with Y=$(ts.properties.Y), Z=$(ts.properties.Z), $(length(ts.AMRs)) EEPs and $(length(ts.properties.masses)) initial stellar mass points.")
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
`p = PARSECLibrary()`, it is callable as
 - `p(Z::Number)` to interpolate the full library to a new metal mass fraction
   (returning a [`PARSECTrackSet`](@ref)), or
 - `p(Z::Number, M::Number)` to interpolate the tracks to a specific metallicity
   and initial stellar mass (returning a [`PARSECTrack`](@ref)).

This type also supports isochrone construction
(see [isochrone](@ref StellarTracks.isochrone(::StellarTracks.PARSEC.PARSECLibrary, ::Number, ::Number))).

# Examples
```jldoctest
julia> p = PARSECLibrary()
Structure of interpolants for PARSEC v1.2S library of stellar tracks. Valid range of metal mass fraction Z is (0.0001, 0.06).

julia> isochrone(p, 10.05, 0.001654) isa NamedTuple
true
```
"""
struct PARSECLibrary{A,B,C} <: AbstractTrackLibrary
    ts::A # Vector of `TrackSet`s
    Z::B  # Vector of Z for each TrackSet
    Y::C  # Vector of Y for each TrackSet
end
# Interpolation to get a TrackSet with metallicity Z
function (ts::PARSECLibrary)(Z::Number)
    error("Not yet implemented.")
end
# Interpolation to get a Track with mass M and metallicity Z
function (ts::PARSECLibrary)(Z::Number, M::Number)
    error("Not yet implemented.")
end
Z(p::PARSECLibrary) = p.Z
Y(p::PARSECLibrary) = p.Y
X(p::PARSECLibrary) = 1 .- p.Y .- p.Z
# MH(p::PARSECLibrary) = PARSEC_MH.(Z(p)) # MH.(Z(p), Y(p))
chemistry(::PARSECLibrary) = PARSECChemistry()
MH(tl::PARSECLibrary) = MH.(chemistry(tl), Z(tl))
Base.eltype(p::PARSECLibrary) = typeof(first(p.Z))
Base.Broadcast.broadcastable(p::PARSECLibrary) = Ref(p)
function Base.show(io::IO, mime::MIME"text/plain", p::PARSECLibrary)
    print(io, "Structure of interpolants for PARSEC v1.2S library of stellar tracks. Valid range of metal mass fraction Z is $(extrema(p.Z)).")
end
function PARSECLibrary(base_dir::AbstractString=datadep"PARSECv1.2S")
    # Load all data into Tables; 1.2s single-threaded, 0.8s multi-threaded
    # ts = [CSV.read(fname, Table) for fname in glob("Z*.gz", base_dir)]
    # Processing into TrackSet structures takes additional time
    set_files = glob("Z*.jld2", base_dir)
    filestems = [splitext(basename(fi))[1] for fi in set_files]
    Z = [parse(Float64, split(split(fi, 'Z')[2], 'Y')[1]) for fi in filestems]
    Y = [parse(Float64, split(fi, 'Y')[2]) for fi in filestems]
    # Sort according to Z; isochrone(p::PARSEC...) depends on this
    idxs = sortperm(Z)
    Z .= Z[idxs]
    Y .= Y[idxs]
    set_files .= set_files[idxs]
    # Make vector of tracksets
    # ts = [PARSECTrackSet(CSV.read(fname, Table), Z[i], Y[i]) for (i, fname) in enumerate(set_files)]
    ts = [PARSECTrackSet(JLD2.load_object(fname), Z[i], Y[i]) for (i, fname) in enumerate(set_files)]
    return PARSECLibrary(ts, Z, Y)
end
"""
    isochrone(p::PARSECLibrary, logAge::Number, Z::Number)
Interpolates properties of the stellar tracks in the library at the requested logarithmic age (`logAge = log10(age [yr])`) and metal mass fraction `Z`. Returns a `NamedTuple` containing the properties listed below:
 - `eep`: Equivalent evolutionary points
 - `m_ini`: Initial stellar masses, in units of solar masses.
 - `logTe`: Base-10 logarithm of the effective temperature [K] of the stellar model.
 - `Mbol`: Bolometric luminosity of the stellar model.
 - `logg`: Surface gravity of the stellar model calculated as `-10.616 + log10(mass) + 4 * logTe - (4.77 - Mbol) / 2.5`.
 - `C_O`: Photospheric C/O ratio (the ZAMS value is used before the TP-AGB).
"""
function isochrone(p::PARSECLibrary, logAge::Number, Z::Number)
    idx = findfirst(Base.Fix1(≈, Z), p.Z) # Will be === nothing if no entry in p.Z is ≈ Z
    # If input Z is represented in base grid, no Z interpolation needed
    if !isnothing(idx) # idx !== nothing
        return isochrone(p.ts[idx], logAge)
    end
    # Check Z is in valid range
    minZ, maxZ = extrema(p.Z)
    if Z < minZ || Z > maxZ
        throw(DomainError(Z, "Requested metallicity Z=$Z is outside the valid range for PARSEC library of $(extrema(p.Z))."))
    end
    # Z is valid, need to interpolate isochrone as a function of Z
    # According to Marigo2017, the interpolations (at least for BCs) in Z or [M/H]
    # are linear, so the BCs overall are computed by a 
    # 3D linear interpolation in logTe x logg x [M/H] space.
    # VandenBerg2014 suggests linear interpolation is sufficient for Z or [M/H] interpolation,
    # although cubic interpolation was used in VandenBerg2012.
    # The pre-computed grid seems more uniform in [M/H] than it is in Z, so I think it might
    # be a good idea to do linear interpolation in [M/H].
    xvec = MH(p) # [M/H] for all TrackSet in p, sorted least to greatest
    x = MH(chemistry(p), Z)  # Requested [M/H]
    # searchsortedfirst returns the index of the first value in xvec greater than or
    # equivalent to x. If x is greater than all values in xvec, return lastindex(xvec) + 1.
    # We have already checked bounds so we know minZ < Z < maxZ
    idx = searchsortedfirst(xvec, x)
    # Evaluate isochrones on either side of intermediate point
    y0 = isochrone(p.ts[idx-1], logAge)
    y1 = isochrone(p.ts[idx], logAge)
    # Get intersection of valid EEPs from each isochrone
    min_eep = max(first(y0.eep), first(y1.eep))
    max_eep = min(last(y0.eep), last(y1.eep))
    # Get indices into y0 and y1 that correspond to the overlapping EEP points
    y0_idxs = Vector{Int}(undef, 0)
    y1_idxs = similar(y0_idxs)
    good_eeps = similar(y0_idxs)
    for (y0_idx, eep) in enumerate(y0.eep)
        y1_idx = searchsortedfirst(y1.eep, eep)
        if y1_idx < lastindex(y1.eep) + 1 # Match found
            push!(y0_idxs, y0_idx)
            push!(y1_idxs, y1_idx)
            push!(good_eeps, eep)
        end
    end
    # Get isochrone keys, removing EEP since that is fixed
    goodkeys = filter(Base.Fix1(!==, :eep), keys(y0))
    # Linear interpolation here
    # result = NamedTuple{goodkeys}((y0[key][y0_idxs] .* (xvec[idx] - x) .+ y1[key][y1_idxs] .* (x - xvec[idx-1])) ./ (xvec[idx] - xvec[idx-1]) for key in goodkeys)
    # Use a function barrier to improve performance since some of the types of the variables
    # aren't known at runtime
    result = NamedTuple{goodkeys}(_interp_kernel(goodkeys, y0, y1, idx, y0_idxs, y1_idxs, x, xvec))
    # Concatenate interpolated result with valid EEP points
    return (eep = good_eeps, result...)
end
_interp_kernel(goodkeys, y0, y1, idx, y0_idxs, y1_idxs, x, xvec) =
    ((y0[key][y0_idxs] .* (xvec[idx] - x) .+ y1[key][y1_idxs] .* (x - xvec[idx-1])) ./ (xvec[idx] - xvec[idx-1]) for key in goodkeys)

#################################################################################

export PARSECLibrary, PARSECChemistry, MH_canon, Z_canon # Unique module exports
export mass, chemistry, X, Y, Z, MH, post_rgb, isochrone # Export generic API methods

#################################################################################

# padova_tracks does some shit to interpolate between the base track and the HB tracks
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

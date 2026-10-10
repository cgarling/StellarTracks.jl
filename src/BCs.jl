"""
    apply_bc(iso, bc::BolometricCorrections.AbstractBCTable)
Given a theoretical isochrone `iso` (e.g., from `isochrone(tracklib, logAge, mh)`) and a bolometric correction table
`bc`, evaluates the bolometric corrections at the effective temperatures and surface gravities of `iso` and applies
them to its bolometric magnitudes, returning a `TypedTables.Table` with the columns of `iso` followed by the magnitude
in each filter of `bc`. Together with [`bc_table`](@ref), this applies several tables (e.g., for different
extinctions) to one isochrone without recomputing it.

```julia
using StellarTracks, BolometricCorrections
p = MISTv1Library(0.0)   # Load MIST library of non-rotating stellar models
m = MISTv1BCGrid("JWST") # Load MIST library of BCs
iso = isochrone(p, 10.0, -1.01)
isos = [apply_bc(iso, bc_table(p, m, -1.01, Av)) for Av in (0.0, 0.5)]
```
"""
function apply_bc(iso, bc::AbstractBCTable)
    # Parsing of `iso` for the arguments of `bc` happens in BolometricCorrections.jl; `bc(iso)` is (filters × stars)
    mags = _parse_Mbol(iso) .- transpose(bc(iso))
    filters = filternames(bc)
    return Table(Table(iso), Table(NamedTuple{filters}(ntuple(k -> view(mags, :, k), Val(length(filters))))))
end

"""
    bc_table(tracklib::AbstractTrackLibrary, bcg::BolometricCorrections.AbstractBCGrid, mh::Number, Av::Number)
Returns the bolometric correction table interpolated from `bcg` for stars of metallicity [M/H] = `mh`, defined for the
chemical mixture of `tracklib`, and V-band extinction `Av`. `mh` is converted to the metal mass fraction `Z` and then
to the [M/H] of the chemical mixture of `bcg`; for grids indexed by [Fe/H] and [α/Fe] (e.g., `MISTv2BCGrid`),
[α/Fe] is taken from `tracklib` (see [`alphaFe`](@ref)). These are the tables used by
`isochrone(tracklib, bcg, logAge, mh, Av)`; see [`apply_bc`](@ref) to apply them to theoretical isochrones.
"""
bc_table(tl::AbstractTrackLibrary, bcg::AbstractBCGrid, mh::Number, Av::Number) = bcg(MH(chemistry(bcg), Z(chemistry(tl), mh)), Av)
function bc_table(tl::AbstractTrackLibrary, bcg::MISTv2BCGrid, mh::Number, Av::Number)
    afe = alphaFe(tl)
    return bcg(FeH(chemistry(bcg), MH(chemistry(bcg), Z(chemistry(tl), mh)), afe), afe, Av)
end

"""
    isochrone(ts::StellarTracks.AbstractTrackSet,
              bc::BolometricCorrections.AbstractBCTable,
              logAge::Number)
Returns an isochrone as a `TypedTables.Table` calculated using the stellar
evolutionary tracks contained in `ts` with bolometric corrections interpolated
from the provided table `bc` at the logarithmic age `logAge`. Column names
can be retrieved with `TypedTables.columnnames`. The result can be converted
to a matrix with `Tables.matrix`.

```julia
using StellarTracks, BolometricCorrections
p = PARSECLibrary()    # Load PARSEC library of stellar models
m = MISTv1BCGrid("JWST") # Load MIST library of BCs
isochrone(p.ts[1], m(MH(p.ts[1]), 0.0), 10.0)
```

```
Table with 35 columns and 1323 rows:
...
```
"""
isochrone(ts::AbstractTrackSet, bc::AbstractBCTable, logAge::Number) =
    apply_bc(isochrone(ts, logAge), bc)
# This generic fallback will work as long as isochrone(tl, logAge, mh) works
"""
    isochrone(tracklib::AbstractTrackLibrary,
              bc::BolometricCorrections.AbstractBCTable,
              logAge::Number,
              mh::Number)
Returns an isochrone as a `TypedTables.Table` calculated using the stellar
evolutionary tracks contained in `tracklib` with bolometric corrections interpolated
from the provided table `bc` at the logarithmic age `logAge` and metallicity
[M/H] = `mh`. Column names can be retrieved with `TypedTables.columnnames`.
The result can be converted to a matrix with `Tables.matrix`.

```julia
using StellarTracks, BolometricCorrections
p = MISTv1Library(0.0)   # Load MIST library of non-rotating stellar models
m = MISTv1BCGrid("JWST") # Load MIST library of BCs
isochrone(p, m(-1.01, 0.0), 10.0, -1.01)
```

```
Table with 36 columns and 1465 rows:
...
```
"""
isochrone(tracklib::AbstractTrackLibrary, bc::AbstractBCTable, logAge::Number, mh::Number) =
    apply_bc(isochrone(tracklib, logAge, mh), bc)

"""
    isochrone(tracklib::AbstractTrackLibrary,
              bcg::AbstractBCGrid, 
              logAge::Number, 
              mh::Number, 
              Av::Number)
Returns an isochrone as a `TypedTables.Table` calculated using the stellar
evolutionary tracks contained in `tracklib` with bolometric corrections interpolated
from the provided bolometric correction grid `bcg` at the logarithmic age 
`logAge` and metallicity [M/H] = `mh`.
Column names can be retrieved with `TypedTables.columnnames`.
The result can be converted to a matrix with `Tables.matrix`.

In general, the stellar evolutionary tracks
and bolometric correction grid may not share the same solar abundance pattern.
In this case, it is common to normalize both to a common total metallicity
quantified by the metal mass fraction `Z` -- for example, this approach is taken
for normalizing the YBC bolometric grid [Chen2019](@citep). This method specifically treats
the input `mh` ([M/H]) as the desired logarithmic abundance for the chemical mixture used
for the stellar evolutionary tracks (`tracklib`), converts this to a metal
mass fraction, then converts that metal mass fraction into the appropriate
logarithmic metallicity [M/H] for the chemical mixture used by the bolometric
correction grid `bcg`. This "corrected" [M/H] value is used when interpolating the
bolometric grid.

If `mh` lies outside the metallicity range of `tracklib`, the stellar tracks are evaluated
at the nearest metallicity in `tracklib` while the bolometric corrections are still evaluated
at `mh`. This is an extrapolation and is similar to what is done by the PARSEC CMD isochrone service.
The idea is that stellar interiors change little at very low metallicity, so combining a [M/H] = -2.2
dex PARSEC isochrone (meaning, interior quantities only) with a [M/H] = -2.5 dex bolometric correction
is a *good enough* approximation of true photometric quantities at [M/H] = -2.5. Whether this is
actually *good enough* depends on your application; if such extrapolation is not acceptable, you
should ensure that `mh` lies within the metallicity range of `tracklib`, for example by checking
that `mh` is between `minimum(MH(tracklib))` and `maximum(MH(tracklib))`. An `mh` outside
that range is an error from `bcg`.

The bolometric correction table is [`bc_table(tl, bcg, mh, Av)`](@ref bc_table), so for grids indexed by [Fe/H] and
[α/Fe] (e.g., `MISTv2BCGrid`), [α/Fe] is taken from `tl` (see [`alphaFe`](@ref)).

```julia
using StellarTracks, BolometricCorrections
p = MISTv1Library(0.0)   # Load MIST library of non-rotating stellar models
m = MISTv1BCGrid("JWST") # Load MIST library of BCs
isochrone(p, m, 10.0, -1.01, 0.0)
```

```
Table with 36 columns and 1465 rows:
...
```
"""
isochrone(tl::AbstractTrackLibrary, bcg::AbstractBCGrid, logAge::Number, mh::Number, Av::Number) =
    isochrone(tl, bc_table(tl, bcg, mh, Av), logAge, clamp(mh, extrema(MH(tl))...))

# This is most useful for constructing large isochrone tables that can then be written out.
# For large grids (3403 isochrones) selecting out a single isochrone takes longer than
# running a new one (9.381 ms vs 1.2 ms). When constructing isochrones in order to make
# partial CMD templates, it is better just to sample them one-by-one in the threaded loop
# rather than trying to pre-generate them all. 
function isochrone(tl::AbstractTrackLibrary,
                   bcg::AbstractBCGrid, logAge::AbstractArray{<:Number}, mh::AbstractArray{<:Number}, Av::Number)
    result = []
    rlock = ReentrantLock()
    # This implementation returns a single Table, with Z and logAge rows 
    Threads.@threads for mhval in mh
        Zval = Z(chemistry(tl), mhval)
        bc = bc_table(tl, bcg, mhval, Av)
        Threads.@threads for la in logAge
            iso = isochrone(tl, bc, la, clamp(mhval, extrema(MH(tl))...))
            # Lock and push into result, adding Z and logAge columns
            @lock rlock push!(result, Table(Table(Z_ini=fill(Zval, length(iso)),
                                                  MH=fill(mhval, length(iso)),
                                                  logAge=fill(la, length(iso))), iso))
        end
    end
    return vcat(result...)
    # This implementation returns a vector of Tables,
    # but you don't necessarily know the order they finished ...
    # zvec = []
    # zlock = ReentrantLock()
    # lavec = []
    # lalock = ReentrantLock()
    # Threads.@threads for Zval in Z
    #     mh = BC.MIST.MH(BC.MIST.MISTv1Chemistry(), Zval)
    #     bc = bcg(mh, Av)
    #     Threads.@threads for la in logAge
    #         iso = isochrone(tl, bc, la, Zval)
    #         @lock rlock push!(result, iso)
    #         @lock zlock push!(zvec, Zval)
    #         @lock lalock push!(lavec, la)
    #     end
    # end
    # # Use zvec and lavec to compute the correct permutation vector to put
    # # result into the correct order. Specifically, ((Z= 0.01, la=(...)), (Z=0.02, la=(...)), ...)
    # # we iterate la more quickly than Z.
    # return result, zvec, lavec
end

# Not sure how to handle the fact that AbstractTrackLibrary and AbstractBCGrid can
# have different dependent variables (Z, Av, α-abundance, etc.). Going to define
# specific call signatures for each type that will be relatively simple to extend.
# For now we are not interpolating against α-abundance and all BC grids have Av
# as an interpolation variable, so the above generic functions are sufficient. 
# Leaving this here as we may need to revisit this when more BC grids are added.

####################################################################################
# Code for specific stellar models / BC grid combinations

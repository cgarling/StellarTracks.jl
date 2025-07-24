"""
    _apply_bc(iso, bc::AbstractBCTable)
Given a calculated isochrone `iso` and an `AbstractBCTable`, parse the isochrone
for `Teff, logg, Mbol`, evaluate BCs as a function of `(Teff, logg)`, and apply
BCs to `Mbol`, returning a joined `TypedTables.Table` with theoretical and
observational quantities.
"""
@inline function _apply_bc(iso, bc::AbstractBCTable)
    iso_Mbol = _parse_Mbol(iso)
    BCs = permutedims(bc(iso)) # Parsing of `iso` for correct args happens in BolometricCorrections.jl
    # Concatenate theoretical quantities with magnitudes
    return Table(Table(iso),
                 Table(Tables.table(iso_Mbol .- BCs; header=filternames(bc))))
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
m = MISTBCGrid("JWST") # Load MIST library of BCs
isochrone(p.ts[1], m(MH(p.ts[1]), 0.0), 10.0)
```

```
Table with 35 columns and 1323 rows:
...
```
"""
isochrone(ts::AbstractTrackSet, bc::AbstractBCTable, logAge::Number) =
    _apply_bc(isochrone(ts, logAge), bc)
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
p = MISTLibrary(0.0)   # Load MIST library of non-rotating stellar models
m = MISTBCGrid("JWST") # Load MIST library of BCs
isochrone(p, m(-1.01, 0.0), 10.0, -1.01)
```

```
Table with 36 columns and 1465 rows:
...
```
"""
isochrone(tracklib::AbstractTrackLibrary, bc::AbstractBCTable, logAge::Number, mh::Number) =
    _apply_bc(isochrone(tracklib, logAge, mh), bc)

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

```julia
using StellarTracks, BolometricCorrections
p = MISTLibrary(0.0)   # Load MIST library of non-rotating stellar models
m = MISTBCGrid("JWST") # Load MIST library of BCs
isochrone(p, m, 10.0, -1.01, 0.0)
```

```
Table with 36 columns and 1465 rows:
...
```
"""
function isochrone(tl::AbstractTrackLibrary,
                   bcg::AbstractBCGrid, logAge::Number, mh::Number, Av::Number)
    # Take stellar track mh, convert to Z, then convert to MH for the BC chemistry
    bc_mh = MH(chemistry(bcg), Z(chemistry(tl), mh))
    return isochrone(tl, bcg(bc_mh, Av), logAge, mh)
end

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
        bc_mh = MH(chemistry(bcg), Zval)
        bc = bcg(bc_mh, Av)
        Threads.@threads for la in logAge
            iso = isochrone(tl, bc, la, mhval)
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
    #     mh = BC.MIST.MH(BC.MIST.MISTChemistry(), Zval)
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
# For now we ware not interpolating against α-abundance and all BC grids have Av
# as an interpolation variable, so the above generic functions are sufficient. 
# Leaving this here as we may need to revisit this when more BC grids are added.

####################################################################################
# Code for specific stellar models / BC grid combinations

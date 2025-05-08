# Isochrones generated from different stellar track libraries may have different keys
# for the same quantities (though we try to avoid this). These functions parse the
# isochrones and return the necessary quantities. 
@inline function _parse_teff(iso)
    iso_keys = keys(iso)
    if :logTe in iso_keys
        return exp10.(iso.logTe)
    elseif :Te in iso_keys
        return iso.Te
    elseif :T in iso
        return iso.T
    else
        throw(ArgumentError("Provided `iso` argument does not contain a recognized effective temperature key, `(:logTe, :Te, :T)`."))
    end 
end

@inline function _parse_logg(iso)
    iso_keys = keys(iso)
    if :logg in iso_keys
        return iso.logg
    elseif :log_g in iso_keys
        return iso.log_g
    else
        throw(ArgumentError("Provided `iso` argument does not contain a recognized surface gravity key, `(:logg, :log_g)`."))
    end
end

@inline function _parse_Mbol(iso)
    iso_keys = keys(iso)
    if :Mbol in iso_keys
        return iso.Mbol
    else
        throw(ArgumentError("Provided `iso` argument does not contain a recognized bolometric magnitude key, `(:Mbol,)`."))
    end
end
"""
    _apply_bc(iso, bc::AbstractBCTable)
Given a calculated isochrone `iso` and an `AbstractBCTable`, parse the isochrone
for `Teff, logg, Mbol`, evaluate BCs as a function of `(Teff, logg)`, and apply
BCs to `Mbol`, returning a joined `TypedTables.Table` with theoretical and
observational quantities.
"""
@inline function _apply_bc(iso, bc::AbstractBCTable)
    iso_teff = _parse_teff(iso)
    iso_logg = _parse_logg(iso)
    iso_Mbol = _parse_Mbol(iso)
    BCs = permutedims(bc(iso_teff, iso_logg))
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

# Examples

```julia
using StellarTracks, BolometricCorrections
p = PARSECLibrary() # Load PARSEC library of stellar models
m = MISTBCGrid()    # Load MIST library of BCs
isochrone(p.ts[1], m(1e-4, 0.0), 10.0)
```

```
Table with 35 columns and 1323 rows:
...
```
"""
isochrone(ts::AbstractTrackSet, bc::AbstractBCTable, logAge::Number) = _apply_bc(isochrone(ts, logAge), bc)
# This generic fallback will work as long as isochrone(tl, logAge, mh) works
isochrone(tl::AbstractTrackLibrary, bc::AbstractBCTable, logAge::Number, mh::Number) =
    _apply_bc(isochrone(tl, logAge, mh), bc)

# Not sure how to handle the fact that AbstractTrackLibrary and AbstractBCGrid can
# have different dependent variables (Z, Av, α-abundance, etc.). Going to define
# specific call signatures for each type that will be relatively simple to extend.

####################################################################################
# Code for PARSEC stellar models

isochrone(tl::PARSECLibrary, bc::AbstractBCTable, logAge::Number, Z::Number) = _apply_bc(isochrone(tl, logAge, Z), bc)
function isochrone(tl::PARSECLibrary, bcg::MISTBCGrid, logAge::Number, Z::Number, Av::Number)
    # MISTBCGrid expects metallicity in [Fe/H] == [M/H], not Z, so convert
    mh = MH(chemistry(bcg), Z)
    return isochrone(tl, bcg(mh, Av), logAge, Z)
end
# This is most useful for constructing large isochrone tables that can then be written out.
# For large grids (3403 isochrones) selecting out a single isochrone takes longer than
# running a new one (9.381 ms vs 1.2 ms). When constructing isochrones in order to make
# partial CMD templates, it is better just to sample them one-by-one in the threaded loop
# rather than trying to pre-generate them all. 
function isochrone(tl::PARSECLibrary, bcg::MISTBCGrid, logAge::AbstractArray{<:Number},
                   Z::AbstractArray{<:Number}, Av::Number)
    result = []
    rlock = ReentrantLock()
    # This implementation returns a single Table, with Z and logAge rows 
    Threads.@threads for Zval in Z
        mh = MH(chemistry(bcg), Zval)
        bc = bcg(mh, Av)
        Threads.@threads for la in logAge
            iso = isochrone(tl, bc, la, Zval)
            # Lock and push into result, adding Z and logAge columns
            @lock rlock push!(result, Table(Table(Z_ini=fill(Zval, length(iso)),
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

####################################################################################
# Code for MIST stellar models

# use generic isochrone(tl::AbstractTrackLibrary, bc::AbstractBCTable, logAge::Number, mh::Number)
isochrone(tl::MISTLibrary, bcg::MISTBCGrid, logAge::Number, mh::Number, Av::Number) =
    isochrone(tl, bcg(mh, Av), logAge, mh)


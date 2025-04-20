module BolometricCorrectionsExt

using StellarTracks: AbstractTrackSet, PARSECLibrary
import StellarTracks: isochrone
using BolometricCorrections: AbstractBCTable, MISTBCGrid, filternames

using TypedTables: Table
import Tables
# BolometricCorrections.MISTBCGrid can be interpolated in metallicity and reddening;
# other grids may contain other

#############################################################

# Isochrones generated from different stellar track libraries may have different keys
# for the same quantities (though we try to avoid this). These functions parse the
# isochrones and return the necessary quantities. 
function _parse_teff(iso)
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

function _parse_logg(iso)
    iso_keys = keys(iso)
    if :logg in iso_keys
        return iso.logg
    else
        throw(ArgumentError("Provided `iso` argument does not contain a recognized surface gravity key, `(:logg,)`."))
    end
end

function _parse_Mbol(iso)
    iso_keys = keys(iso)
    if :Mbol in iso_keys
        return iso.Mbol
    else
        throw(ArgumentError("Provided `iso` argument does not contain a recognized bolometric magnitude key, `(:Mbol,)`."))
    end
end

"""
    isochrone([::Type{TypedTables.Table},]
              ts::StellarTracks.AbstractTrackSet,
              bc::BolometricCorrections.AbstractBCTable,
              logAge::Number)
Returns an isochrone as a matrix with dimensions `(nstars, nfilters)` calculated using the stellar evolutionary tracks contained in `ts` with bolometric corrections interpolated from the provided table `bc` at the logarithmic age `logAge`. The three-argument version that returns a `Table` has a roughly fixed runtime overhead cost of 3--5 μs to perform the type conversion.

# Examples
Default call signature returning a matrix.
```julia
using StellarTracks, BolometricCorrections
p = PARSECLibrary() # Load PARSEC library of stellar models
m = MISTBCGrid()    # Load MIST library of BCs
isochrone(p.ts[1], m(1e-4, 0.0), 10.0)
```

```
1323×29 Matrix{Float64}
...
```

Call signature returning a `TypedTables.Table`.
```julia
using TypedTables: Table
isochrone(Table, p.ts[1], m(1e-4, 0.0), 10.0)
```

```
Table with 29 columns and 1323 rows:
...
```
"""
function isochrone(ts::AbstractTrackSet, bc::AbstractBCTable, logAge::Number)
    iso = isochrone(ts, logAge)
    iso_teff = _parse_teff(iso)
    iso_logg = _parse_logg(iso)
    iso_Mbol = _parse_Mbol(iso)
    BCs = permutedims(bc(iso_teff, iso_logg))
    return iso_Mbol .- BCs
end
function isochrone(::Type{Table}, ts::AbstractTrackSet, bc::AbstractBCTable, logAge::Number)
    filters = filternames(bc)
    # Mostly fixed ~3--5 μs runtime cost
    return Table(Tables.table(isochrone(ts, bc, logAge); header=filters))
end

# Not sure how to handle the fact that AbstractTrackLibrary and AbstractBCGrid can
# have different dependent variables (Z, A_V, α-abundance, etc.). For now
# going to set up some specific call signatures that be relatively simple to extend.
function isochrone(tl::PARSECLibrary, bc::AbstractBCTable, logAge::Number, Z::Number)
    iso = isochrone(tl, logAge, Z)
    iso_teff = _parse_teff(iso)
    iso_logg = _parse_logg(iso)
    iso_Mbol = _parse_Mbol(iso)
    BCs = permutedims(bc(iso_teff, iso_logg))
    return iso_Mbol .- BCs
end
function isochrone(::Type{Table}, tl::PARSECLibrary, bc::AbstractBCTable,
                   logAge::Number, Z::Number)
    filters = filternames(bc)
    # Mostly fixed ~3--5 μs runtime cost
    return Table(Tables.table(isochrone(tl, bc, logAge, Z); header=filters))
end
isochrone(tl::PARSECLibrary, bcg::MISTBCGrid, logAge::Number, Z::Number, A_V::Number) =
    isochrone(tl, bcg(Z, A_V), logAge, Z)
isochrone(::Type{Table}, tl::PARSECLibrary, bcg::MISTBCGrid,
          logAge::Number, Z::Number, A_V::Number) =
    isochrone(Table, tl, bcg(Z, A_V), logAge, Z)

end # module

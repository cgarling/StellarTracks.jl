using StellarTracks, BolometricCorrections
using TypedTables: Table

p = PARSECLibrary()
m = MISTBCGrid("JWST")

# Currently can't interpolate a new TrackSet from a PARSECLibrary
# instance, so we just take the first TrackSet from p.ts
iso = isochrone(p.ts[1], m(1e-4, 0.0), 10.0)
@test iso isa Matrix{Float64}
iso_table = isochrone(Table, p.ts[1], m(1e-4, 0.0), 10.0)
@test iso_table isa Table

# Test passing full stellar track library with single BC table
iso1 = isochrone(p, m(1e-4, 0.0), 10.0, 1e-4)
@test iso1 isa Matrix{Float64}
# Test passing full stellar track library with full BC grid
iso2 = isochrone(p, m, 10.0, 1e-4, 0.0)
@test iso2 isa Matrix{Float64}
# These should both be equivalent
@test iso1 == iso2
# Test table interface
t1 = isochrone(Table, p, m(1e-4, 0.0), 10.0, 1e-4)
@test t1 isa Table
t2 = isochrone(Table, p, m, 10.0, 1e-4, 0.0)
@test t2 isa Table
# These should both be equivalent
@test t1 == t2

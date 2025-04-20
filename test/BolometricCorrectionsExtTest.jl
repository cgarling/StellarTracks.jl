using StellarTracks: isochrone, PARSECLibrary, Z
using BolometricCorrections: MISTBCGrid
using TypedTables: Table
using Test

p = PARSECLibrary()
m = MISTBCGrid("JWST")

# Currently can't interpolate a new TrackSet from a PARSECLibrary
# instance, so we just take the first TrackSet from p.ts
iso = isochrone(p.ts[1], m(1e-4, 0.0), 10.0)
@test iso isa Table

# Test passing full stellar track library with single BC table
iso1 = isochrone(p, m(1e-4, 0.0), 10.0, 1e-4)
@test iso1 isa Table
# Test passing full stellar track library with full BC grid
iso2 = isochrone(p, m, 10.0, 1e-4, 0.0)
@test iso2 isa Table
# These should both be equivalent
@test iso1 == iso2

# Test that we can interpolate across the range of valid Z
for ZZ in range(extrema(Z(p))...; length=10)
    @test isochrone(p, m, 10.0, ZZ, 0.0) isa Table
end

# Test that we can interpolate across the range of valid A_V
# for ZZ in range(extrema(Z(p))...; length=10)
#     @test_nowarn isochrone(p, m, 10.0, ZZ, 0.0)
# end

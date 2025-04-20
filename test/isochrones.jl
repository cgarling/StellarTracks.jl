using StellarTracks

p = PARSECLibrary()
# Test isochrone(ts::StellarTracks.AbstractTrackSet, logAge::Number)
iso1 = isochrone(p.ts[1], 10.0)
@test iso1 isa NamedTuple
iso2 = isochrone(p, 10.0, Z(first(p.ts)))
# Could do a regression test here
@test iso2 isa NamedTuple
@test keys(iso1) == keys(iso2)
# iso1 and iso2 should be identical, since no Z interpolation is occuring
for key in keys(iso1)
    @test getproperty(iso1, key) == getproperty(iso2, key)
end
# Test that we can interpolate across the range of valid Z
for ZZ in range(extrema(Z(p))...; length=10)
    @test_nowarn isochrone(p, 10.0, ZZ)
end
# Make sure error is thrown for out of Z range
@test_throws DomainError isochrone(p, 10.0, 1e-11)
@test_throws DomainError isochrone(p, 10.0, 1.0)


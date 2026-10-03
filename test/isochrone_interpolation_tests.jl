using StellarTracks
using Test

struct PairingTrackSet <: StellarTracks.AbstractTrackSet
    data::NamedTuple
end
StellarTracks.isochrone(ts::PairingTrackSet, ::Number) = ts.data

struct PairingLibrary <: StellarTracks.AbstractTrackLibrary
    ts::Vector{PairingTrackSet}
end
StellarTracks.MH(::PairingLibrary) = [-1.0, 0.0]

# Blend only identical EEPs, not the next available point across gaps in either isochrone.
# Check interpolated properties, native-grid results, and disjoint EEP sets.
@testset "Exact EEP pairing across metallicities" begin
    lower = (eep=[1, 2, 3, 4, 5, 6, 8], m_ini=Float64[1, 2, 3, 4, 5, 6, 8], logTe=Float64[11, 12, 13, 14, 15, 16, 18])
    upper = (eep=[2, 4, 5, 7], m_ini=Float64[20, 40, 50, 70], logTe=Float64[120, 140, 150, 170])
    lib = PairingLibrary([PairingTrackSet(lower), PairingTrackSet(upper)])
    iso = isochrone(lib, 10.0, -0.75)
    @test iso.eep == [2, 4, 5]
    @test iso.m_ini == [6.5, 13.0, 16.25]
    @test iso.logTe == [39.0, 45.5, 48.75]
    @test isochrone(lib, 10.0, -1.0) === lower
    @test isochrone(lib, 10.0, 0.0) === upper

    matched = PairingLibrary([PairingTrackSet(lower), PairingTrackSet(lower)])
    @test isochrone(matched, 10.0, -0.5) == lower
    disjoint = PairingLibrary([PairingTrackSet(lower), PairingTrackSet((eep=[7], m_ini=[70.0], logTe=[170.0]))])
    @test all(isempty, values(isochrone(disjoint, 10.0, -0.5)))
end

using HexEarth
import GeoInterface as GI
using Test

@testset "HexEarth.jl" begin
    @testset "Cell Basics" begin
        ll = LatLon((35.9099875, -79.0752876)...)
        cell10 = Cell(ll, 10)
        for res in 11:15
            cell = Cell(ll, res)
            ll2 = GI.centroid(cell)
            @test cell10 != cell
            @test ll.lat ≈ ll2[1] atol=1e-3
            @test ll.lon ≈ ll2[2] atol=1e-3
            @test HexEarth.resolution(cell) == res
            @test HexEarth.is_cell(cell)
        end
    end
end

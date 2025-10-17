using HexEarth
import GeoInterface as GI
using Test

@testset "HexEarth.jl" begin
    @testset "Geo" begin
        a = (0.0, 0.0)
        for angle in 0:359, dist in (10, 100, 1000, 1_000_000)
            b = HexEarth.destination(a, angle, dist)
            @test HexEarth.haversine(a, b) ≈ dist atol=1e-3
            @test HexEarth.bearing(a, b) ≈ angle atol=1e-3
        end
        for a in [rand(2) for _ in 1:1000]
            ex = HexEarth.point2extent(a, 1000)
            x2 = HexEarth.destination(a, 0, 1000)[2]
            @test ex.Y[2] ≈ x2 atol=1e-3
        end
    end

    @testset "Cell Basics" begin
        cell0 = Cell((-79.0752876, 35.9099875), 0)
        c = GI.centroid(cell0)
        for res in 1:15
            cell = Cell(c, res)
            c2 = GI.centroid(cell)
            @test is_cell(cell)
            @test resolution(cell) == res
            @test cell0 != cell
            @test c[1] ≈ c2[1] atol=1e-3
            @test c[2] ≈ c2[2] atol=1e-3
        end
    end
    @testset "Cell Indexing" begin
        a = Cell((0,0))
        @test a[0,0] == a
        @test a[0,0,0] == a
        @test a[1,1][-1,-1] == a
        @test a[1,1,1][-1,-1,-1] == a
        b = a[1, 1]
        @test b != a
        @test resolution(b) == resolution(a)
        @test b[-1, -1] == a
    end
    @testset "cells" begin
        xy = (0,0)
        @test length(cells(xy)) == 1

        # Point
        @test length(cells(GI.Point((0,0)))) == 1

        # MultiPoint
        @test length(cells(GI.MultiPoint([(0,0), (1,1), (2,2)]))) == 3

        # Line
        a = Cell((0,0))
        b = a[1,1]
        line = GI.Line([GI.centroid(a), GI.centroid(b)])
        @test length(cells(line)) == 2

        # LineString
        linestring = GI.LineString([GI.centroid(a), GI.centroid(b), GI.centroid(a[1,0])])
        @test length(cells(linestring)) == 3

        # PolyGon
        poly = GI.Polygon([GI.LineString([GI.centroid(a[1,1]), GI.centroid(a[1,0]), GI.centroid(a[0,0]), GI.centroid(a[0,1]), GI.centroid(a[1,1])])])
        @test length(cells(poly)) == 4

    end
end

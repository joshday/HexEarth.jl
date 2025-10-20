module HexEarth

using H3
using Statistics
using StyledStrings: @styled_str

import GeoInterface as GI
import GeoFormatTypes as GFT
import GeometryOps as GO
import Extents

import H3.API as API
import H3.Lib: LatLng


export
    Cell, Vertex, DirectedEdge,
    cells, datacells, ring, grid_distance, resolution, is_cell, is_vertex, is_directed_edge, is_pentagon

#-----------------------------------------------------------------------------# Notes
# In this package:
#  - All coordinates are assumed to be in EPSG:4326 (WGS84 lat/lon in degrees)
#  - latitude/longitude units == degrees
#  - distance units == meters
#  - bearing units == degrees clockwise from North

#-----------------------------------------------------------------------------# includes
# Essentially all things not hexagon-related
include("Geo.jl")
import .Geo: destination, haversine, bearing, point2extent

# # Hexagon types and functions.  Generated from https://www.redblobgames.com/grids/hexagons/
# include("Hexagons.jl")
# import .Hexagons: ring


#-----------------------------------------------------------------------------# H3.Lib.LatLng
GI.isgeometry(::LatLng) = true
GI.geomtrait(::LatLng) = GI.PointTrait()
GI.coordinates(::GI.PointTrait, o::LatLng) = (rad2deg(o.lng), rad2deg(o.lat))
GI.getcoord(::GI.PointTrait, o::LatLng, i::Integer) = GI.coordinates(o)[i]
GI.ncoord(::GI.PointTrait, ::LatLng) = 2

#-----------------------------------------------------------------------------# H3.API wrappers
"""
    resolution(::UInt64)
    resolution(::Cell)
    resolution(::Vertex)
    resolution(::DirectedEdge)

Return the resolution 0-15 (coarse to fine) of the given H3 index.
"""
resolution(index::UInt64) = API.getResolution(index)

is_cell(index::UInt64) = API.isValidCell(index)
is_directed_edge(index::UInt64) = API.isValidDirectedEdge(index)
is_vertex(index::UInt64) = API.isValidVertex(index)
is_pentagon(index::UInt64) = API.isPentagon(index)

check(f, x::UInt64) = f(x) ? x : throw(ArgumentError("Invalid H3 index: $x does not satisfy condition $f"))

#-----------------------------------------------------------------------------# split_index
"""
    split_index(x::UInt64) -> NamedTuple

Split a 64-bit H3 index into its fields per the spec:

- `mode` (4 bits)
- `resolution` (4 bits, 0–15)
- `basecell` (7 bits, 0–121)
- `digits` (15 three-bit digits; unused digits are 7)

References:
  - H3 bit layout & cell index fields. See docs.
"""
@inline function split_index(x::UInt64)
    mode       = UInt8((x >> 59) & 0xF)    # 4 bits
    resolution = UInt8((x >> 52) & 0xF)    # 4 bits
    basecell   = UInt8((x >> 45) & 0x7F)   # 7 bits

    # Digit i (1-based) lives at bits (45 - 3i) .. (47 - 3i)
    digits = ntuple(i -> UInt8((x >> (3*(15 - i))) & 0x7), 15)

    return (mode=mode, resolution=resolution, basecell=basecell, digits=digits)
end

#-----------------------------------------------------------------------------# H3IndexType
"One of the H3 index types: Cell, Vertex, or DirectedEdge."
abstract type H3IndexType end

GI.isgeometry(::H3IndexType) = true
GI.geomtrait(::H3IndexType) = GI.UnknownTrait()
GI.crs(::H3IndexType) = GFT.EPSG(4326)
GI.ncoord(::GI.PolygonTrait, o::H3IndexType) = 2

is_cell(o::H3IndexType) = is_cell(o.index)
is_directed_edge(o::H3IndexType) = is_directed_edge(o.index)
is_vertex(o::H3IndexType) = is_vertex(o.index)
resolution(o::H3IndexType) = resolution(o.index)
is_pentagon(o::H3IndexType) = is_pentagon(o.index)

#-----------------------------------------------------------------------------# Cell
"""
    Cell(index::UInt64)
    Cell(lonlat::Tuple{<:Real, <:Real}, res::Integer=10)

A hexagonal cell in the H3 grid.
"""
struct Cell <: H3IndexType
    index::UInt64
    Cell(x::UInt64) = new(check(is_cell, x))
end
Cell(lonlat, res=10) = Cell(API.latLngToCell(API.LatLng(deg2rad(lonlat[2]), deg2rad(lonlat[1])), Int(res)))

GI.geomtrait(::Cell) = GI.PolygonTrait()
GI.centroid(::GI.PolygonTrait, o::Cell) = (ll = API.cellToLatLng(o.index); (rad2deg(ll.lng), rad2deg(ll.lat)))
GI.coordinates(::GI.PolygonTrait, o::Cell) = (out = GI.coordinates.(API.cellToBoundary(o.index)); return [out..., out[1]])
GI.nhole(::GI.PolygonTrait, o::Cell) = 0
GI.ngeom(::GI.PolygonTrait, o::Cell) = 1
GI.getgeom(::GI.PolygonTrait, o::Cell, i::Integer) = GI.LineString(GI.coordinates(o))
GI.area(::GI.PolygonTrait, o::Cell) = area(o)  # in m²

function GI.extent(::GI.PolygonTrait, o::Cell)
    coords = GI.coordinates(o)
    lons = getindex.(coords, 1)
    lats = getindex.(coords, 2)
    Extents.Extent(X=(minimum(lons), maximum(lons)), Y=(minimum(lats), maximum(lats)))
end

function Base.show(io::IO, o::Cell)
    shape = is_pentagon(o) ? styled"{bright_red:⬠}" : styled"{bright_green:⬡}"
    ll = styled"{bright_black:$(GI.centroid(o))}"
    print(io, styled"$shape {bright_cyan:$(typeof(o))} {bright_magenta:$(resolution(o))} {bright_black:$ll}")
end

# Check if the cell crosses the 180° longitude line
function crosses_lon180(o::Cell)
    lons = [x[1] for x in GI.coordinates(o)]
    for x in lons, y in lons
        abs(x - y) > 180 && return true
    end
    return false
end

vertices(o::Cell) = Vertex.(API.cellToVertexes(o.index))
const vertexes = vertices

area(o::Cell) = API.cellAreaM2(o.index)  # in m²

parent(o::Cell) = Cell(API.cellToParent(o.index, resolution(o) - 1))
children(o::Cell) = Cell.(API.cellToChildren(o.index, resolution(o) + 1))


grid_distance(a::Cell, b::Cell) = API.gridDistance(a.index, b.index)

"""
    grid_path_cells(a::Cell, b::Cell)

Return a minimal-length contiguous path of cells from `a` to `b`.
"""
grid_path_cells(a::Cell, b::Cell) = Cell.(API.gridPathCells(a.index, b.index))

"""
    grid_disk(o::Cell, k::Integer)

Return all cells within `k` steps of origin cell `o`, including `o` itself.
"""
grid_disk(o::Cell, k::Integer) = Cell.(API.gridDisk(o.index, k))


"""
    grid_ring_unsafe(o::Cell, k::Integer)

Return all cells exactly `k` steps from origin cell `o`.
The order of the returned cells is not guaranteed.

Per h3 documentation, "This function may fail if pentagonal distortion is encountered."  We don't know what that means.
"""
grid_ring_unsafe(o::Cell, k::Integer) = Cell.(API.gridRingUnsafe(o.index, k))

haversine(a::Cell, b::Cell) = haversine(GI.centroid(a), GI.centroid(b))

destination(a::Cell, azimuth°, m) = Cell(destination(GI.centroid(a), azimuth°, m), resolution(a))

#-----------------------------------------------------------------------------# Cell Indexing
const neighbor_indices = ((0,1,0), (1,1,0), (1,0,0), (1,0,1), (0,0,1), (0,1,1))


Base.getindex(o::Cell, i::Integer) = i < 7 ? getindex(o, neighbor_indices[i]...) : throw(BoundsError(o, i))
Base.iterate(o::Cell, i=1) = i > 6 ? nothing : (o[i], i + 1)
Base.IteratorSize(::Type{Cell}) = Base.HasLength()
Base.length(::Cell) = 6
Base.eltype(::Type{Cell}) = Cell

Base.getindex(o::Cell, i::Integer, j::Integer) = GridIJ(o)[i, j]

Base.getindex(o::Cell, i::Integer, j::Integer, k::Integer) = GridIJK(o)[i, j, k]

#-----------------------------------------------------------------------------# cells (geom-to-Vector{Cell})
"""
    cells(geometry, res = 10)
    cells(polygon, res = 10; containment = nothing)
    cells(raster, res = 10; dropmissing = true, dropempty = true, containment = :overlap)

Return a `Vector{Cell}` covering the given geometry at the specified H3 resolution (default 10).

- (Multi)polygons allow specifying a `containment` mode:
    - `nothing`: (default) use libh3's default polygonToCells behavior (we think this is equivalent to `:center`)
    - `:center`: include cells whose center point is within the polygon
    - `:full`: include cells fully contained within the polygon
    - `:overlap`: include cells that overlap the polygon at all
    - `:overlap_bbox`: include cells that overlap the polygon's bounding box
- Rasters (requires Rasters.jl to be loaded)
"""
cells(geom, res::Integer = 10; kw...) = cells(GI.trait(geom), geom, res; kw...)

cells(::Nothing, x::AbstractVector{<:Tuple{<:Real, <:Real}}, res::Integer) = unique!([Cell(x, res) for x in x])

cells(trt::GI.PointTrait, geom, res::Integer) = [Cell(GI.coordinates(trt, geom), res)]

function cells(trait::GI.MultiPointTrait, geom, res::Integer)
    unique!(Cell.(GI.coordinates(trait, geom), res))
end

ϵ::Float64 = 1e-6

function cells(trait::GI.LineTrait, geom, res::Integer; shortest_path = true)
    coords = GI.coordinates(trait, geom)
    a = Cell(coords[1], res)
    b = Cell(coords[2], res)
    out = grid_path_cells(a, b)  # grid_path_cells gives approximate line
    # Now expand to include neighboring cells that the line crosses
    shortest_path || for cell in out, candidate in cell
        candidate in out && continue
        GO.disjoint(geom, candidate) && continue
        push!(out, candidate)
    end
    out
end

function cells(trait::GI.LineStringTrait, geom, res::Integer; shortest_path=true)
    coords = GI.coordinates(trait, geom)
    out = Cell[]
    @views for (a,b) in zip(coords[1:end-1], coords[2:end])
        line = GI.Line([a, b])
        union!(out, cells(line, res; shortest_path=shortest_path))
    end
    return out
end

# Get H3.Lib.GeoPolygon from GeoInterface polygon
function h3polygon(geom)
    GI.trait(geom) == GI.PolygonTrait() || throw(ArgumentError("Expected Polygon geometry"))
    verts = map(GI.coordinates(geom)) do ring
        map(ring) do coord
            API.LatLng(deg2rad(coord[2]), deg2rad(coord[1]))
        end
    end
    GC.@preserve verts begin
        geo_loops = H3.Lib.GeoLoop.(length.(verts), pointer.(verts))
    end
    GC.@preserve geo_loops begin
        n = length(geo_loops)
        return H3.Lib.GeoPolygon(geo_loops[1], n - 1, n > 1 ? pointer(geo_loops, 2) : C_NULL)
    end
end

function cells(trait::GI.PolygonTrait, geom, res::Integer; containment = nothing)
    isnothing(containment) || containment in (:center, :full, :overlap, :overlap_bbox) ||
        throw(ArgumentError("Invalid containment mode.  Expected one of `(nothing, :center, :full, :overlap, :overlap_bbox)`.  Found: $containment."))
    geo_polygon = h3polygon(geom)
    flag = containment == :center ? H3.Lib.CONTAINMENT_CENTER :
        containment == :full ? H3.Lib.CONTAINMENT_FULL :
        containment == :overlap ? H3.Lib.CONTAINMENT_OVERLAPPING :
        containment == :overlap_bbox ? H3.Lib.CONTAINMENT_OVERLAPPING_BBOX :
        H3.Lib.CONTAINMENT_INVALID  # containment = nothing --> this is unused
    GC.@preserve geo_polygon begin
        max_n = Ref{Int64}()
        ret::API.H3Error = isnothing(containment) ?
            H3.Lib.maxPolygonToCellsSize(Ref(geo_polygon), res, 0, max_n) :
            H3.Lib.maxPolygonToCellsSizeExperimental(Ref(geo_polygon), res, flag, max_n)
        API._check_h3error(ret, nothing)
        out = zeros(UInt64, max_n[])
        ret2::API.H3Error = isnothing(containment) ?
            H3.Lib.polygonToCells(Ref(geo_polygon), res, 0, out) :
            H3.Lib.polygonToCellsExperimental(Ref(geo_polygon), res, flag, max_n[], out)
        API._check_h3error(ret2, out)
        out = Cell.(filter!(!iszero, unique!(out)))
    end
end

function cells(trait::GI.MultiPolygonTrait, geom, res::Integer; kw...)
    reduce(union, cells.(GI.getpolygon(geom), res; kw...))
end

function cells(ex::Extents.Extent, res::Integer; kw...)
    x1, x2 = ex.X
    y1, y2 = ex.Y
    ls = GI.LineString([(x1, y1), (x1, y2), (x2, y2), (x2, y1), (x1, y1)])
    cells(GI.Polygon([ls]), res; kw...)
end

#-----------------------------------------------------------------------------# GridIJ
"""
    GridIJ(origin::Cell)

Create a 2D grid of cells indexed by (i, j) coordinates relative to an origin cell.  Axes are 120° apart.

- Unlike libh3, (0, 0) will always correspond to the origin cell.
- Pentagon distortion may result in some cells having multiple (i, j) coordinates.
- See also `GridIJK`.
"""
struct GridIJ
    origin::Cell
    ij::API.CoordIJ
end
GridIJ(o::Cell) = GridIJ(o, API.cellToLocalIj(o.index, o.index))
Base.show(io::IO, o::GridIJ) = print(io, styled"GridIJ - origin: $(o.origin)")

function Base.getindex(grid::GridIJ, i::Integer, j::Integer)
    Cell(API.localIjToCell(grid.origin.index, API.CoordIJ(i + grid.ij.i, j + grid.ij.j)))
end

function Base.getindex(grid::GridIJ, o::Cell)
    ij = API.cellToLocalIj(grid.origin.index, o.index)
    return (ij.i - grid.ij.i, ij.j - grid.ij.j)
end

#-----------------------------------------------------------------------------# GridIJK
"""
    GridIJK(origin::Cell)

Create a 2D grid of cells indexed by (i, j, k) coordinates relative to an origin cell.

- Unlike libh3, (0, 0, 0) will always correspond to the origin cell.
- Note there are multiple valid (i, j, k) coordinates for a given cell.
- See also `GridIJ`.
"""
struct GridIJK
    origin::Cell
    ijk::API.CoordIJK
end
GridIJK(o::Cell) = GridIJK(o, API.cellToLocalIjk(o.index, o.index))
Base.show(io::IO, o::GridIJK) = print(io, styled"GridIJK - origin: $(o.origin)")
function Base.getindex(grid::GridIJK, i::Integer, j::Integer, k::Integer)
    coord = API.CoordIJK(i + grid.ijk.i, j + grid.ijk.j, k + grid.ijk.k)
    res = API.localIjkToCell(grid.origin.index, coord)
    res isa Union{H3.API.H3ErrorCode, H3.API.H3Error} ?
        error(H3.API.describeH3Error(res)) :
        Cell(res)
end

#-----------------------------------------------------------------------------# Vertex
"""
    Vertex(index::UInt64)

Represents a vertex (point on a `Cell` boundary) in the H3 grid.
"""
struct Vertex <: H3IndexType
    index::UInt64
    Vertex(x::UInt64) = new(check(is_vertex, x))
end
GI.geomtrait(::Vertex) = GI.PointTrait()
GI.coordinates(::GI.PointTrait, o::Vertex) = (ll = API.vertexToLatLng(o.index); (rad2deg(ll.lng), rad2deg(ll.lat)))
GI.getcoord(::GI.PointTrait, o::Vertex, i::Integer) = GI.coordinates(o)[i]

function Base.show(io::IO, o::Vertex)
    ll = styled"{bright_black:$(GI.coordinates(o))}"
    print(io, styled"{bright_green:•} {bright_cyan:$(typeof(o))} {bright_magenta:$(resolution(o))} {bright_black:$(repr(o.index))} $ll")
end

#-----------------------------------------------------------------------------# DirectedEdge
"""
    DirectedEdge(a::Cell, b::Cell)

Represents a directional relationship between *adjacent* cells `a` and `b`.
"""
struct DirectedEdge <: H3IndexType
    index::UInt64
    DirectedEdge(x::UInt64) = new(check(is_directed_edge, x))
end
function DirectedEdge(a::Cell, b::Cell)
    idx = API.cellsToDirectedEdge(a.index, b.index)
    idx isa H3.API.H3Error ?
        error(H3.API.describeH3Error(idx)) :
        DirectedEdge(idx)
end

GI.geomtrait(::DirectedEdge) = GI.LineTrait()
GI.coordinates(::GI.LineTrait, o::DirectedEdge) = GI.centroid.(cells(o))
GI.getcoord(::GI.LineTrait, o::DirectedEdge, i::Integer) = GI.coordinates(o)[i]

function Base.show(io::IO, o::DirectedEdge)
    ll = styled"{bright_black:$(GI.coordinates(o))}"
    print(io, styled"{bright_green:→} {bright_cyan:$(typeof(o))} {bright_magenta:$(resolution(o))} {bright_black:$(repr(o.index))} $ll")
end

function cells(o::DirectedEdge)
    a, b = API.directedEdgeToCells(o.index)
    [Cell(a), Cell(b)]
end

Base.length(o::DirectedEdge) = 2
Base.getindex(o::DirectedEdge, i::Integer) = cells(o)[i]
Base.iterate(o::DirectedEdge, i=1) = iterate(cells(o), i)
Base.IteratorSize(::Type{DirectedEdge}) = Base.HasLength()
Base.eltype(::Type{DirectedEdge}) = Cell

end # module

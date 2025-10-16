module HexEarth

using H3
using H3_jll
using StyledStrings: @styled_str

import GeoInterface as GI
import GeoFormatTypes as GFT
import Rasters as R
import H3.API as API
import Extents

export
    LatLon, Cell, Vertex,
    cells, geocode, resolution, is_cell, is_vertex, is_directed_edge, is_pentagon

#-----------------------------------------------------------------------------# Notes
# In this package:
#  - latitude/longitude units == degrees
#  - distance units == meters
#  - bearing units == degrees clockwise from North

#-----------------------------------------------------------------------------# LatLon
"""
    LatLon(lat::T, lon::T) where {T<:Real}
    LatLon(; lat, lon)
    LatLon(::Tuple{<:Real, <:Real})
    LatLon(::AbstractVector{<:Real})

A simple struct to represent a latitude/longitude pair in *degrees* (EPSG:4326).
"""
struct LatLon{T}
    lat::T
    lon::T
end
LatLon(o::API.LatLng) = LatLon(rad2deg(o.lat), rad2deg(o.lng))
LatLon(; lat, lon) = LatLon(lat, lon)
LatLon(x::AbstractVector{<:Real}) = length(x) == 2 ? LatLon(x[2], x[1]) : throw(ArgumentError("Expected 2-element vector, got length $(length(x))"))
LatLon(x::Tuple{<:Real, <:Real}) = LatLon(x[2], x[1])

Base.NamedTuple((; lon, lat)::LatLon) = (x=lon, y=lat)
Base.show(io::IO, o::LatLon) = print(io, NamedTuple(o))

H3.Lib.LatLng(o::LatLon) = H3.Lib.LatLng(deg2rad(o.lat), deg2rad(o.lon))

GI.isgeometry(::LatLon) = true
GI.geomtrait(::LatLon) = GI.PointTrait()
GI.coordinates(::GI.PointTrait, o::LatLon) = (o.lon, o.lat)
GI.getcoord(::GI.PointTrait, o::LatLon, i::Integer) = GI.coordinates(o)[i]
GI.ncoord(::GI.PointTrait, ::LatLon) = 2

"Approximate radius of the earth in meters (WGS84)."
const R = 6_371_000

"""
    haversine(a::LatLon, b::LatLon)
Calculate the great-circle distance (meters) between two LatLon points using the Haversine formula.
"""
function haversine(a::LatLon{T}, b::LatLon{T}) where {T <: Real}
    x = sind((b.lat - a.lat) / 2) ^ 2 + cosd(a.lat) * cosd(b.lat) * sind((b.lon - a.lon) / 2) ^ 2
    return 2R * asin(min(sqrt(x), one(x)))
end

"""
    destination(a::LatLon, bearing°, dist_m)

Find destination point given starting point (LatLon), bearing (clockwise from North), and distance (m)
"""
function destination(a::LatLon, bearing°, m)
    ϕ1, λ1 = a.lat, a.lon
    δ = rad2deg(m / R)
    ϕ2 = asind(sind(ϕ1) * cosd(δ) + cosd(ϕ1) * sind(δ) * cosd(bearing°))
    λ2 = λ1 + atand(sind(bearing°) * sind(δ) * cosd(ϕ1), cosd(δ) - sind(ϕ1) * sind(ϕ2))
    LatLon(ϕ2, λ2)
end


#-----------------------------------------------------------------------------# extent
"""
    Extents.extent(o::LatLon; n=1000, s=1000, e=1000, w=1000)
    Extents.extent(o::LatLon, meters)

Create a bounding box (Extents.Extent) around a LatLon point `o`, extending `n`, `s`, `e`, and `w` meters in each direction .

"""
function Extents.extent(o::LatLon; n=1000, s=1000, e=1000, w=1000)
    N = destination(o, 0, n).lat
    E = destination(o, 90, e).lon
    S = destination(o, 180, s).lat
    W = destination(o, 270, w).lon
    return Extents.Extent(X=(W, E), Y=(S, N))
end
Extents.extent(o::LatLon, meters) = Extents.extent(o; n=meters, s=meters, e=meters, w=meters)

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
struct Cell <: H3IndexType
    index::UInt64
    Cell(x::UInt64) = new(check(is_cell, x))
end
Cell(o::LatLon, resolution=10) = Cell(API.latLngToCell(API.LatLng(deg2rad(o.lat), deg2rad(o.lon)), Int(resolution)))
Cell(x::AbstractString) = Cell(geocode(x))

GI.geomtrait(::Cell) = GI.PolygonTrait()
GI.centroid(::GI.PolygonTrait, o::Cell) = LatLon(API.cellToLatLng(o.index))
GI.coordinates(::GI.PolygonTrait, o::Cell) = (out = LatLon.(API.cellToBoundary(o.index)); return [out..., out[1]])
GI.nhole(::GI.PolygonTrait, o::Cell) = 0
GI.ngeom(::GI.PolygonTrait, o::Cell) = 1
GI.getgeom(::GI.PolygonTrait, o::Cell, i::Integer) = GI.LineString(GI.coordinates(o))
GI.area(::GI.PolygonTrait, o::Cell) = area_m2(o)  # in m²

function Base.show(io::IO, o::Cell)
    shape = is_pentagon(o) ? styled"{bright_red:⬠}" : styled"{bright_green:⬡}"
    ll = styled"{bright_black:$(GI.centroid(o))}"
    print(io, styled"$shape {bright_cyan:$(typeof(o))} {bright_magenta:$(resolution(o))} {bright_black:$(repr(o.index))} $ll")
end

vertices(o::Cell) = Vertex.(API.cellToVertexes(o.index))
const vertexes = vertices

area_m2(o::Cell) = API.cellAreaM2(o.index)  # in m²

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

destination(a::Cell, bearing°, m) = Cell(destination(GI.centroid(a), bearing°, m), resolution(a))


Base.getindex(o::Cell, i::Integer) = Vertex(API.cellToVertex(o.index, i))

Base.getindex(o::Cell, i::Integer, j::Integer) = GridIJ(o)[i, j]

Base.getindex(o::Cell, i::Integer, j::Integer, k::Integer) = GridIJK(o)[i, j, k]


#-----------------------------------------------------------------------------# cells (geom-to-Vector{Cell})
"""
    cells(geometry, res = 10)

Return a `Vector{Cell}` covering the given geometry at the specified H3 resolution (default 10).
"""
cells(geom, res::Integer = 10; kw...) = cells(GI.trait(geom), geom, res; kw...)

cells(trait::GI.PointTrait, geom, res::Integer) = [Cell(LatLon(GI.coordinates(trait, geom)...), res)]

function cells(trait::GI.MultiPointTrait, geom, res::Integer)
    unique!(Cell.(LatLon.(GI.coordinates(trait, geom), res)))
end

function cells(trait::GI.LineStringTrait, geom, res::Integer)
    coords = GI.coordinates(trait, geom)
    out = [Cell(LatLon(reverse(coords[1])...), res)]
    for (lon, lat) in @view coords[2:end]
        c = Cell(LatLon(lat, lon), res)
        path = grid_path_cells(out[end], c)
        append!(out, path[2:end])
    end
    unique!(out)
end

function cells(trait::GI.PolygonTrait, geom, res::Integer)
    verts = [H3.Lib.LatLng.(LatLon.(ring)) for ring in GI.coordinates(trait, geom)]
    GC.@preserve verts begin
        geo_loops = H3.Lib.GeoLoop.(length.(verts), pointer.(verts))
    end
    GC.@preserve geo_loops begin
        n = length(geo_loops)
        geo_polygon = H3.Lib.GeoPolygon(geo_loops[1], n - 1, n > 1 ? pointer(geo_loops, 2) : C_NULL)
    end
    GC.@preserve geo_polygon begin
        max_n = Ref{Int64}()
        ret::API.H3Error = H3.Lib.maxPolygonToCellsSize(Ref(geo_polygon), res, 0, max_n)
        API._check_h3error(ret, nothing)
        out = zeros(UInt64, max_n[])
        ret2::API.H3Error = H3.Lib.polygonToCells(Ref(geo_polygon), res, 0, out)
        API._check_h3error(ret2, out)
        out = Cell.(filter!(!iszero, unique!(out)))
    end
    # Hack around the fact that the libh3 implementation sometimes misses cells on the edges.
    for ring in GI.coordinates(geom)
        union!(out, cells(GI.LineString(ring), res))
    end
    return out
end

function cells(trait::GI.MultiPolygonTrait, geom, res::Integer)
    reduce(union, cells.(GI.getpolygon(geom), res))
end

function cells(ex::Extents.Extent, res::Integer)
    sw = LatLon(ex.Y[1], ex.X[1])
    ne = LatLon(ex.Y[2], ex.X[2])
    nw = LatLon(ex.Y[2], ex.X[1])
    se = LatLon(ex.Y[1], ex.X[2])
    cells(GI.Polygon([GI.LineString([sw, nw, ne, se, sw])]), res)
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
"""
struct GridIJK
    origin::Cell
    ijk::API.CoordIJK
end
GridIJK(o::Cell) = GridIJK(o, API.cellToLocalIjk(o.index, o.index))
Base.show(io::IO, o::GridIJK) = print(io, styled"GridIJK - origin: $(o.origin)")
function Base.getindex(grid::GridIJK, i::Integer, j::Integer, k::Integer)
    Cell(API.localIjkToCell(grid.origin.index, API.CoordIJK(i + grid.ijk.i, j + grid.ijk.j, k + grid.ijk.k)))
end

#-----------------------------------------------------------------------------# Vertex
struct Vertex <: H3IndexType
    index::UInt64
    Vertex(x::UInt64) = new(check(is_vertex, x))
end
GI.geomtrait(::Vertex) = GI.PointTrait()
GI.coordinates(::GI.PointTrait, o::Vertex) = LatLon(API.vertexToLatLng(o.index))
GI.getcoord(::GI.PointTrait, o::Vertex, i::Integer) = GI.coordinates(o)[i]

function Base.show(io::IO, o::Vertex)
    ll = styled"{bright_black:$(GI.coordinates(o))}"
    print(io, styled"{bright_cyan:$(typeof(o))} {bright_magenta:$(resolution(o))} {bright_black:$(repr(o.index))} $ll")
end

LatLon(o::Vertex) = LatLon(API.vertexToLatLng(o.index))


#-----------------------------------------------------------------------------# Makie interop


end # module

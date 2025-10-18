```@setup geom
using GLMakie, HexEarth


```
# Geometry to Cell

The main function to create a collection of cells that approximate a geometry is [`cells`](@ref):

```julia
cells(geometry, res::Int)::Vector{Cell}
```

## Point

```@example geom
pt = (0, 0)

cell = Cell(pt)

p = lines(cell)
scatter!(pt, color=0)
p
```

## LineString

- Notice that the lines are not entirely covered by hexagons.  The `cells` function will create a
shortest path between the two points and does not check intermediate points on the line.


```@example geom
import GeoInterface as GI

ls = GI.LineString([(0,0), (1,1), (-1,2)])

x = cells(ls, 6)

p = lines(x)
lines!(ls.geom)
p
```

## (Multi)Polygon

```@example geom
using GeoJSON

# Get GeoJSON polygon of North Carolina
url = "https://nominatim.openstreetmap.org/search?state=north_carolina&format=geojson&polygon_geojson=1"
file = download(url)
obj = GeoJSON.read(file)

x = cells(obj.geometry[1], 6)

p = lines(x)
lines!(obj.geometry)
p
```

## Rasters

- For things that have data associated with points, such as rasters, `cells` returns `DataCells`.
- `DataCells{T}` is a wrapper around a `Dict{Cell, T}`.
- For rasters, multiple points could end up in the same cell, so we use `DataCells{Vector{T}}` where `T` is the eltype of the raster.
- To reduce vector values to a single number, use `HexEarth.get_values(::DataCell, fun)`.

```@example geom
using Rasters, RasterDataSources, ArchGDAL

(; elev) = getraster(WorldClim{Elevation})

r = Raster(elev)
madagascar = view(r, X(43.25 .. 50.48), Y(-25.61 .. -12.04))

x = cells(madagascar, 4)

poly(x; color = HexEarth.get_values(x, maximum))
```

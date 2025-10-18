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
```

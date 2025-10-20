```@setup geom
using CairoMakie, HexEarth
```

# Geometry to Cells

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

- By default, the `cells` function will create a shortest path between the cells that contain the two points of each line.  You can set `shortest_path=false` to include all hexagons that overlap the lines.


```@example geom
import GeoInterface as GI

ls = GI.LineString([(0,0), (1,1), (-1,2)])

x = cells(ls, 6)
x2 = cells(ls, 6; shortest_path=false)

fig = Figure()
ax = Axis(fig[1,1], title="shortest_path = true")
ax2 = Axis(fig[1,2], title="shortest_path = false")
lines!(ax, x)
lines!(ax, ls.geom)
lines!(ax2, x2)
lines!(ax2, ls.geom)
fig
```

## (Multi)Polygon

- There are several algorithms available for filling in cells based on a polygon, based on keyword argument `containment`.
- `containment = nothing` (the default) will use the `H3.Lib.polygonToCells` function.
- Other containment modes (`:center`, `:full`, `:overlap`, or `:overlap_bbox`) will use `H3.Lib.polygonToCellsExperimental`.
- As far as we can tell `containment = :center` uses the same algorithm as `containment = nothing`.

```@example geom
using GeoJSON

# Get GeoJSON polygon of North Carolina
url = "https://nominatim.openstreetmap.org/search?state=north_carolina&format=geojson&polygon_geojson=1"
file = download(url)
obj = GeoJSON.read(file)

# Create figure
function make_axis(i, j, res, containment)
    ax = Axis(fig[i, j], title="containment = $containment")
    x = cells(obj.geometry[1], res; containment)
    lines!(ax, x)
    lines!(ax, obj.geometry)
end
fig = Figure()
make_axis(1, 1, 4, nothing)
make_axis(1, 2, 4, :center)
make_axis(1, 3, 4, :full)
make_axis(2, 1, 4, :overlap)
make_axis(2, 2, 4, :overlap_bbox)
fig
```

## Extents

- Uses the same `containment` options as for (Multi)polygons.

```@example geom
ex = GI.extent(obj)

# for plotting
x1, x2 = ex.X
y1, y2 = ex.Y
linestring = GI.LineString([(x1,y1), (x1,y2), (x2,y2), (x2, y1), (x1,y1)])

# Create figure
function make_axis(i, j, res, containment)
    ax = Axis(fig[i, j], title="containment = $containment")
    x = cells(ex, res; containment)
    # lines!(ax, obj.geometry)
    lines!(ax, linestring)
    lines!(ax, x)
end
fig = Figure()
make_axis(1, 1, 4, :center)
make_axis(1, 2, 4, :full)
make_axis(2, 1, 4, :overlap)
make_axis(2, 2, 4, :overlap_bbox)
fig
```

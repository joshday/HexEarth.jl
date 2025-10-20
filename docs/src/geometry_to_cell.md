```@setup geom
using CairoMakie, HexEarth
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

- The `cells` function has more control over how cells are placed for polygon inputs (see below).  You can hack together a polygon from a linestring by adding some noise to the points in reverse order and closing the loop.

```@example geom
ls = GI.LineString([(0.0, 0.0), (1.0, 1.0), (-1.0, 2.0)])

# Create linestring for polygon
out = [
    # Original points
    GI.coordinates(ls)...,
    # Reverse points with added noise
    map(xy -> xy .+ .0000001, reverse(GI.coordinates(ls)))...,
    # Close the loop
    GI.coordinates(ls)[1]
]
shape = GI.Polygon([out])

# include cells that have any overlapping part with the line
x = cells(shape, 6, containment = :overlap)

p = lines(x)
lines!(ls)
p
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

## Rasters

- For rasters, `cells` returns `Dict{Cell, Vector{T}}` where `T` is the eltype of the raster.

```@example geom
using Rasters, RasterDataSources, ArchGDAL

(; elev) = getraster(WorldClim{Elevation})

r = Raster(elev)

nc_ext = GI.extent(obj.geometry[1])

r_nc = mask(r, with=obj.geometry[1])

nc_elev = view(r_nc, nc_ext)

x = cells(nc_elev, 4)

poly(collect(keys(x)); color=maximum.(values(x)))
```

### Resolution Issues

- At higher `Cell` resolutions, you may end up with holes (cells for which no raster pixel is contained).

```@example geom
x = cells(nc_elev, 5)

poly(collect(keys(x)); color=maximum.(values(x)))
```

- `cells` accepts a "mask" as a first argument which will:
  1.  Create the cells based on the geometry.
  2.  Populate the cells' data based on the raster.
- Cells that would otherwise have been empty will now use the nearest raster pixel.

```@example geom
x = cells(obj.geometry[1], nc_elev, 5)

poly(collect(keys(x)); color=maximum.(values(x)))
```

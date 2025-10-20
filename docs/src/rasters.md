# Rasters to Cells

- For rasters, `cells` returns `Dict{Cell, Vector{T}}` where `T` is the eltype of the raster.

```@example geom
using HexEarth, GeoMakie, CairoMakie, Rasters, RasterDataSources, ArchGDAL
import GeoInterface as GI
import GeometryOps as GO

r = Raster(getraster(WorldClim{Elevation}).elev)

dict = cells(r, 2)

# remove hexagons overlapping 180° longitude to avoid plotting
# horizontal "stripes" across the map
dict2 = filter(kv -> !HexEarth.crosses_lon180(kv[1]), dict)

poly(collect(keys(dict2)), color=maximum.(values(dict2)), axis=(;type = GeoAxis))
```

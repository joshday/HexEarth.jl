# HexEarth.jl

HexEarth provides tools built on top of [H3.jl](https://github.com/wookay/H3.jl) for working with [Uber's H3 grid system](https://h3geo.org).

!!! tip "What does HexEarth do?"
    The primary purpose of this package is to join disparate vector and raster datasets to a unified hexagonal grid.

## Quick Start

- All coordinates are assumed to be in EPSG:4326 (WGS84 lon/lat in degrees)
- Hexagon cell resolutions are 0-15 (coarse-to-fine).  Resolution 15 corresponds to cells with areas smaller than 1 meter.
- All distances are assumed to be (and returned in) meters.

```@repl
using HexEarth
import GeoInterface as GI

new_york = (-75.0, 43.0)

cell = Cell(new_york, 10)

HexEarth.area(cell)  # meters ^ 2

GI.centroid(cell)
```

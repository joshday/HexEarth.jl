module HexEarthMakieExt

using HexEarth
using Makie

import GeoInterface as GI

GI.@enable_makie Makie HexEarth.H3IndexType
GI.@enable_makie Makie HexEarth.LatLon

end

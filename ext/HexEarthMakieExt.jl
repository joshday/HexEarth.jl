module HexEarthMakieExt

using Makie, HexEarth
import GeoInterface as GI

GI.@enable_makie Makie H3IndexType
GI.@enable_makie Makie LatLon

end

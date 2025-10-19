module HexEarthMakieExt

import Makie
import HexEarth: H3IndexType, DataCells
import GeoInterface as GI

GI.@enable_makie Makie H3IndexType

end

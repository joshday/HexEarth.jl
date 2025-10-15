module HexEarthGMTExt

import GMT
import HexEarth: geocode, LatLon

geocode(x::AbstractString) = LatLon(GMT.geocoder(x).data...)

end

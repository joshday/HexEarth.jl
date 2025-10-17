module Geo

import GeoInterface as GI
import Extents

export haversine, destination, bearing, point2extent

# distances in meters
# angles in degrees
# coordinate order is (lon, lat)

#-----------------------------------------------------------------------------# R (earth radius)
"Approximate radius of the earth in meters (WGS84)."
const R = 6_371_000

#-----------------------------------------------------------------------------# haversine
"""
    haversine(a, b)

Calculate the great-circle distance (meters) between two (lon, lat) points using the Haversine formula.
`a` and `b` are any iterables of length 2 representing (lon, lat).
"""
function haversine((a_lon, a_lat), (b_lon, b_lat))
    x = sind((b_lat - a_lat) / 2) ^ 2 + cosd(a_lat) * cosd(b_lat) * sind((b_lon - a_lon) / 2) ^ 2
    return 2R * asin(min(sqrt(x), one(x)))
end

#-----------------------------------------------------------------------------# destination
"""
    destination(x, azimuth°, dist_m)

Find destination point given starting point `x` (longitude/latitude in degrees), `azimuth` (clockwise from North), and `dist_m` (meters).
"""
function destination((lon, lat), azimuth°, m)
    δ = rad2deg(m / R)
    lat2 = asind(sind(lat) * cosd(δ) + cosd(lat) * sind(δ) * cosd(azimuth°))
    lon2 = lon + atand(sind(azimuth°) * sind(δ) * cosd(lat), cosd(δ) - sind(lat) * sind(lat2))
    (lon2, lat2)
end

#-----------------------------------------------------------------------------# bearing
function bearing((a_lon, a_lat), (b_lon, b_lat))              # assume degrees in, degrees out
    Δlon = b_lon - a_lon

    # atan2(y, x) with degree trig
    y = cosd(b_lat) * sind(Δlon)
    x = cosd(a_lat)*sind(b_lat) - sind(a_lat)*cosd(b_lat)*cosd(Δlon)

    θ = atand(y, x)                      # [-180, 180]
    return mod(θ + 360, 360)             # [0, 360)
end

#-----------------------------------------------------------------------------# point2extent
point2extent(x, dist) = point2extent(x; n=dist, s=dist, e=dist, w=dist)

function point2extent(x; n=1000, s=1000, e=1000, w=1000)
    N = destination(x, 0, n)[2]
    E = destination(x, 90, e)[1]
    S = destination(x, 180, s)[2]
    W = destination(x, 270, w)[1]
    return Extents.Extent(X=(W, E), Y=(S, N))
end

end

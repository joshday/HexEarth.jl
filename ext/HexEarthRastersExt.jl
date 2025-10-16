module HexEarthRastersExt

import Rasters as R
import HexEarth as H
import GeoInterface as GI

# cells(r::R.AbstractRaster, res::Integer=10; dropmissing::Bool=true) = cells(GI.RasterTrait(), r, res; dropmissing)

function H.cells(trait::GI.RasterTrait, geom, res::Integer; dropmissing::Bool = true)
    x = R.lookup(r, R.X)
    y = R.lookup(r, R.Y)
    out = dropmissing ? Dict{H.Cell, Vector{Base.nonmissingtype(T)}}() : Dict{H.Cell, Vector{T}}()

    for i in 1:size(r, 1), j in 1:size(r, 2)
        val = r[i, j]
        if dropmissing && ismissing(val)
            # Skip missing values
        else
            lat = x[i]
            lon = y[j]
            cell = H.Cell(H.LatLon(lat, lon), resolution)
            v = get!(out, cell, T[])
            push!(v, val)
        end
    end

    return out
end

end

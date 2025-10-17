module HexEarthRastersExt

import Rasters
import HexEarth: cells, LatLon, Cell


function cells(r::Rasters.AbstractRaster{T, 2}, res::Integer = 10; dropmissing::Bool = true) where {T}
    x = Rasters.lookup(r, Rasters.X)
    y = Rasters.lookup(r, Rasters.Y)
    out = dropmissing ? Dict{Cell, Vector{Base.nonmissingtype(T)}}() : Dict{Cell, Vector{T}}()

    for i in eachindex(x), j in eachindex(y)
        val = r[i, j]
        if dropmissing && ismissing(val)
            # Skip missing values
        else
            lon = x[i]
            lat = y[j]
            cell = Cell((lat, lon), res)
            v = get!(out, cell, T[])
            push!(v, val)
        end
    end

    return out
end

end

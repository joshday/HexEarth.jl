module HexEarthRastersExt

import Rasters: AbstractRaster, X, Y, At
import HexEarth: cells, Cell, DataCells
import GeoInterface as GI


function cells(r::AbstractRaster{T, 2}, res::Integer = 10; dropmissing::Bool = true) where {T}
    out = dropmissing ? Dict{Cell, Vector{Base.nonmissingtype(T)}}() : Dict{Cell, Vector{T}}()

    for lon in r.dims[1], lat in r.dims[2]
        val = r[X(At(lon)), Y(At(lat))]
        if !(dropmissing && ismissing(val))
            cell = Cell((lon, lat), res)
            v = get!(out, cell, T[])
            push!(v, val)
        end
    end
    return DataCells(out)
end

end

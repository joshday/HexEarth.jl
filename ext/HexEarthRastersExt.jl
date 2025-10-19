module HexEarthRastersExt

import Rasters: AbstractRaster, X, Y, At, Near
import HexEarth: cells, Cell
import GeoInterface as GI


function cells(r::AbstractRaster{T, 2}, res::Integer; dropmissing::Bool = true) where {T}
    out = dropmissing ? Dict{Cell, Vector{Base.nonmissingtype(T)}}() : Dict{Cell, Vector{T}}()

    for lon in r.dims[1], lat in r.dims[2]
        val = r[X(At(lon)), Y(At(lat))]
        if !(dropmissing && ismissing(val))
            cell = Cell((lon, lat), res)
            v = get!(out, cell, T[])
            push!(v, val)
        end
    end
    return out
end

function cells(mask, r::AbstractRaster{T, 2}, res::Integer; dropmissing::Bool = true) where {T}
    x = cells(mask, res)
    out = Dict{Cell, Vector{T}}()
    for cell in x
        lon, lat = GI.centroid(cell)
        val = r[X(Near(lon)), Y(Near(lat))]
        v = get!(out, cell, T[])
        if ismissing(val) && dropmissing
            delete!(out, cell)
        else
            push!(v, val)
        end
    end
    return out
end

end

module HexEarthRastersExt

import Rasters: AbstractRaster, X, Y, At, Near, dims
import HexEarth: cells, bind, Cell
import GeoInterface as GI

function cells(r::AbstractRaster{T, 2}, res::Integer = 10; dropmissing = true, containment = :overlap) where {T}
    # Step 1) Initialize output hexagons to cover the extent of the raster
    all_cells = cells(GI.extent(r), res; containment)
    S = dropmissing ? Base.nonmissingtype(T) : T
    out = Dict(c => S[] for c in all_cells)

    # Step 2) Iterate through the raster and populate hexagons
    for lon in r.dims[1], lat in r.dims[2]
        val = r[X(At(lon)), Y(At(lat))]
        cell = Cell((lon, lat), res)
        v = get!(out, cell, T[])  # get! is required depending on containment mode
        (!ismissing(val) || !dropmissing) && push!(v, val)
    end

    # Step 3) For cells without data, find the closest point from the raster
    for (k, v) in out
        if isempty(v)
            lon, lat = GI.centroid(k)
            val = r[X(Near(lon)), Y(Near(lat))]
            (!ismissing(val) || !dropmissing) && push!(v, val)
        end
    end

    # Step 4) drop empty vectors (only happens for missing values)
    dropmissing && filter!(kv -> !isempty(kv[2]), out)
    return out
end

end

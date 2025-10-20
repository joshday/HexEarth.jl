# Indexing

- Cells are uniquely identified by a `UInt64` *index*.  The only global indexing system is through these values.
  - You can create cells from its index or from coordinates that are contained in the cell: `Cell(index)` or `Cell((lon, lat), res=10)`.
- However, H3 supports several other indexing systems that are valid within a single face of the icosahedron.
- HexEarth uses several different `getindex` methods to expose these indexing systems.

!!! important "Index Ordering"
    Index ordering is not consistent.  I.e. there are no directional guarantees as to where `cell[1]` is
    relative to `cell`.

```@repl inspect
cell[1]  # spiral indexing (1-6) for immediate neighbors

cell[0, 1]  # https://h3geo.org/docs/core-library/coordsystems#local-ij-coordinates

cell[0, 1, 0]  # https://h3geo.org/docs/core-library/coordsystems#ijk-coordinates
```

### Spiral Indexing

```@example inspect
fig = Figure()
ax = Axis(fig[1,1])
poly!(ax, cell)
for i in 1:6
    lines!(ax, cell[i])
    text!(ax, GI.centroid(cell[i]), text=string(i))
end
fig
```

### IJ Indexing

!!! note
    For IJ/IJK indexing, HexEarth uses the conventions that (0,0) and (0,0,0) refer to the origin cell.  For whatever reason, this is not the case in the original h3 library.

- Reference: [https://h3geo.org/docs/core-library/coordsystems/#local-ij-coordinates](https://h3geo.org/docs/core-library/coordsystems/#local-ij-coordinates)

```@example inspect
fig = Figure()
ax = Axis(fig[1,1])
poly!(ax, cell)
for i in [(0,1), (0,-1), (-1,-1), (1,1), (1,0), (-1,0)]
    lines!(ax, cell[i...])
    text!(ax, GI.centroid(cell[i...]), text=string(i))
end
fig
```

### IJK Indexing

- Reference: [https://h3geo.org/docs/core-library/coordsystems#ijk-coordinates](https://h3geo.org/docs/core-library/coordsystems#ijk-coordinates)
- IJK addresses are not unique unless "normalized"
- Note: normalized IJK coordinates have at most two non-zero coordinates.

![](https://h3geo.org/images/ijkp.png)

```@example inspect
fig = Figure()
ax = Axis(fig[1,1])
poly!(ax, cell)
for i in HexEarth.neighbor_indices
    @info i
    lines!(ax, cell[i...])
    text!(ax, GI.centroid(cell[i...]), text=string(i))
end
fig
```

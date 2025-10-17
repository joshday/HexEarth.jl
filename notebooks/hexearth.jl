### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 560f27e0-aac5-11f0-8a6d-497eb8fb50cb
begin 
	using Pkg, Revise
	Pkg.activate(joinpath(@__DIR__))
	
	using HexEarth
	import GeoInterface as GI 
	using GLMakie, GeoMakie, Extents, Rasters, RasterDataSources, ArchGDAL
end

# ╔═╡ ea82065b-fa89-4e20-97ee-042988e10eaf
md"# Indexing"

# ╔═╡ 850c51c3-9bf8-4cc5-a2db-c2459e3c4d1d
o = Cell((40.7128, -74.006))

# ╔═╡ 606d6716-fc91-4353-b560-632217c3d66e
o[1]

# ╔═╡ 151368a5-9de2-4e3d-bdc4-120637a9682f
let 
	fig = Figure()
	ax = Axis(fig[1, 1])
	scatter!(ax, o)
	lines!(ax, o)
	scatter!(ax, [o[i] for i in 1:6])
	fig
end

# ╔═╡ 71de79a8-7dac-436c-bd77-2baf2a4f1339
o[1]

# ╔═╡ b7c18070-209f-4f4d-af37-8ae4271e7142
md"""
# Geometry to Cells
"""

# ╔═╡ 36c1fe57-7caf-42f7-bb23-de5a8ba28136
md"## Point"

# ╔═╡ c83218a0-50fe-4f52-ba97-60132974021c
point = LatLon(-74.006, 40.7128)

# ╔═╡ 126171cf-9aab-4c47-933b-b2fea37f1f07
point_cells = cells(point)

# ╔═╡ f5907824-e8f8-4959-b92f-8db18246ed31
let 
	fig = Figure()
	ax = Axis(fig[1,1])
	poly!(point_cells)
	scatter!(point; color=:black)
	fig
end

# ╔═╡ ac25ae60-3cf1-4396-b36f-682e93fe98a9
md"## MultiPoint"

# ╔═╡ 751998f4-5ee7-491c-9f7d-0dc1225f5c76
multipoint = [point, HexEarth.destination(point, 45, 1000)]

# ╔═╡ e09bf92d-0ab9-40c5-916e-8ef63392d425
multipoint_cells = cells(GI.MultiPoint(multipoint))

# ╔═╡ 223901ff-a4dd-4b6d-9daa-9fdb28a41e50
let 
	fig = Figure()
	ax = Axis(fig[1,1])
	poly!(ax, multipoint_cells)
	scatter!(ax, multipoint, color=:black)
	fig
end

# ╔═╡ 4ff6fda0-1b55-478c-b94f-721ee82d478a
md"## LineString"

# ╔═╡ 66de4d2a-3542-41d4-a3c6-579f2f3e60c5
linestring = [HexEarth.destination(point, i, 500) for i in (0:60:360) .+ 15]

# ╔═╡ 04c49090-09b2-4054-8bb7-2c60336abc51
linestring_cells = cells(GI.LineString(linestring), 11)

# ╔═╡ c9fa611a-bb9f-4ba7-91eb-845518ea5988
let 
	fig = Figure()
	ax = Axis(fig[1,1])
	poly!(ax, linestring_cells)
	lines!(ax, linestring_cells)
	scatter!(ax, linestring, color=:black)
	lines!(ax, linestring, color=:black)
	fig
end

# ╔═╡ eeef9192-71b2-4442-9d08-95949815c326
md"## Polygon"

# ╔═╡ 7579c136-6cb9-4432-a502-9ee8be762773
polygon_cells = cells(GI.Polygon([linestring]), 11)

# ╔═╡ 92b07948-4b9d-4edf-8670-2ad405947089
let 
	fig = Figure()
	ax = Axis(fig[1,1])
	poly!(ax, polygon_cells)
	lines!(ax, polygon_cells)
	lines!(ax, linestring, color=:black)
	fig
end

# ╔═╡ db180a74-3d30-4934-a534-88d256222b1b
md"## MultiPolygon"

# ╔═╡ a2a88eaa-9300-4678-a88d-16b8dc4931d4
multipolygon = [[linestring], [HexEarth.destination.(linestring, 17, 1500)]]

# ╔═╡ 255310d0-ef0b-4f74-b97e-0f3396045aea
multipolygon_cells = cells(GI.MultiPolygon(multipolygon), 11)

# ╔═╡ ebca2f56-9dd6-40fa-b91c-8f1366938ce9
let 
	fig = Figure()
	ax = Axis(fig[1,1])
	poly!(ax, multipolygon_cells)
	lines!(ax, multipolygon_cells)
	poly!(ax, GI.MultiPolygon(multipolygon), alpha=.5)
	fig
end

# ╔═╡ db0df0d5-1256-48e6-8048-374d44bb7469
md"## Extent"

# ╔═╡ b590af1b-fd93-4c8d-a864-77d5119b2ccc
ex = Extents.Extent(X=(0.0, 0.05), Y=(0.0, 0.05))

# ╔═╡ f26dcc93-7c96-41de-a3eb-38e9a618313a
ex_cells = cells(ex)

# ╔═╡ 1298c152-d52f-4830-9045-9567e301387e
let 
	fig = Figure()
	ax = Axis(fig[1,1])
	poly!(ax, ex_cells)
	lines!(ax, ex_cells)
	ex_poly = [[
		(ex.X[1], ex.Y[1]),
		(ex.X[1], ex.Y[2]),
		(ex.X[2], ex.Y[2]),
		(ex.X[2], ex.Y[1]),
		(ex.X[1], ex.Y[1]),
	]]
	poly!(ax, GI.Polygon(ex_poly), color=:black, alpha=.3)
	fig
end

# ╔═╡ a421b4f7-1d58-44bc-8585-9c3af3070315
md"## Raster"

# ╔═╡ d1eaedd6-f2c4-45ea-bf69-2415bc180f38
r = Raster(getraster(WorldClim{Climate}, :wind; month=1))

# ╔═╡ 14a39fcd-5c02-4e86-88e3-23d3e0dfaea3
r_cells = cells(r, 2)

# ╔═╡ 1a94de8e-ceb7-474d-b179-860b9986bd7e
let 
	fig = Figure()
	ax = Axis(fig[1,1])
	ax2 = Axis(fig[2,1])
	plot!(ax, r)
	poly!(ax2, collect(keys(r_cells)))
	fig
end

# ╔═╡ Cell order:
# ╟─560f27e0-aac5-11f0-8a6d-497eb8fb50cb
# ╟─ea82065b-fa89-4e20-97ee-042988e10eaf
# ╠═850c51c3-9bf8-4cc5-a2db-c2459e3c4d1d
# ╠═606d6716-fc91-4353-b560-632217c3d66e
# ╠═151368a5-9de2-4e3d-bdc4-120637a9682f
# ╠═71de79a8-7dac-436c-bd77-2baf2a4f1339
# ╟─b7c18070-209f-4f4d-af37-8ae4271e7142
# ╟─36c1fe57-7caf-42f7-bb23-de5a8ba28136
# ╟─c83218a0-50fe-4f52-ba97-60132974021c
# ╟─126171cf-9aab-4c47-933b-b2fea37f1f07
# ╟─f5907824-e8f8-4959-b92f-8db18246ed31
# ╟─ac25ae60-3cf1-4396-b36f-682e93fe98a9
# ╟─751998f4-5ee7-491c-9f7d-0dc1225f5c76
# ╟─e09bf92d-0ab9-40c5-916e-8ef63392d425
# ╟─223901ff-a4dd-4b6d-9daa-9fdb28a41e50
# ╟─4ff6fda0-1b55-478c-b94f-721ee82d478a
# ╠═66de4d2a-3542-41d4-a3c6-579f2f3e60c5
# ╟─04c49090-09b2-4054-8bb7-2c60336abc51
# ╟─c9fa611a-bb9f-4ba7-91eb-845518ea5988
# ╟─eeef9192-71b2-4442-9d08-95949815c326
# ╠═7579c136-6cb9-4432-a502-9ee8be762773
# ╟─92b07948-4b9d-4edf-8670-2ad405947089
# ╟─db180a74-3d30-4934-a534-88d256222b1b
# ╟─a2a88eaa-9300-4678-a88d-16b8dc4931d4
# ╟─255310d0-ef0b-4f74-b97e-0f3396045aea
# ╟─ebca2f56-9dd6-40fa-b91c-8f1366938ce9
# ╟─db0df0d5-1256-48e6-8048-374d44bb7469
# ╟─b590af1b-fd93-4c8d-a864-77d5119b2ccc
# ╟─f26dcc93-7c96-41de-a3eb-38e9a618313a
# ╟─1298c152-d52f-4830-9045-9567e301387e
# ╟─a421b4f7-1d58-44bc-8585-9c3af3070315
# ╠═d1eaedd6-f2c4-45ea-bf69-2415bc180f38
# ╠═14a39fcd-5c02-4e86-88e3-23d3e0dfaea3
# ╠═1a94de8e-ceb7-474d-b179-860b9986bd7e

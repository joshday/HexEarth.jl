using Documenter
using HexEarth
using Makie
Makie.set_theme!(theme_dark())

get!(ENV, "RASTERDATASOURCES_PATH", mktempdir())

makedocs(
    sitename = "HexEarth",
    format = Documenter.HTML(),
    modules = [HexEarth],
    pages = [
        "Home" => "index.md",
        "Geometry to Cells" => "geometry_to_cell.md",
        "Rasters to Cells" => "rasters.md",
        "Inspection and Indexing" => "inspection.md",
        "API" => "api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/joshday/HexEarth.jl",
)

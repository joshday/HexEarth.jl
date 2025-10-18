using Documenter
using HexEarth
using Makie
Makie.set_theme!(theme_dark())

makedocs(
    sitename = "HexEarth",
    format = Documenter.HTML(),
    modules = [HexEarth]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/joshday/HexEarth.jl",
)

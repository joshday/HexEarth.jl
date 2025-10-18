# Claude-generated module, based on https://www.redblobgames.com/grids/hexagons
module Hexagons
#-----------------------------------------------------------------------------# Coordinate Types
"""
    Cube(q, r, s)

Cube coordinates for hexagonal grids. Must satisfy the constraint: q + r + s = 0.
This coordinate system is most versatile for hexagon grid algorithms.
"""
struct Cube{T <: Integer}
    q::T
    r::T
    s::T

    function Cube(q::T, r::T, s::T) where T <: Integer
        q + r + s == 0 || error("Cube coordinates must satisfy q + r + s = 0, got $(q + r + s)")
        new{T}(q, r, s)
    end
end

Cube(q, r) = Cube(q, r, -q-r)

"""
    Axial(q, r)

Axial coordinates for hexagonal grids. Equivalent to cube coordinates with implicit s = -q-r.
More compact storage than cube coordinates while maintaining most of the benefits.
"""
struct Axial{T <: Integer}
    q::T
    r::T
end

"""
    OffsetCoord(col, row)

Offset coordinates for hexagonal grids. Used with offset types: OddR, EvenR, OddQ, EvenQ.
"""
struct OffsetCoord{T<:Integer}
    col::T
    row::T
end

# Offset coordinate types
abstract type OffsetType end
struct OddR <: OffsetType end
struct EvenR <: OffsetType end
struct OddQ <: OffsetType end
struct EvenQ <: OffsetType end

#-----------------------------------------------------------------------------# Conversions
Base.convert(::Type{Cube}, a::Axial) = Cube(a.q, a.r, -a.q-a.r)
Base.convert(::Type{Axial}, c::Cube) = Axial(c.q, c.r)

# Offset coordinate conversions
function offset_to_axial(::Type{OddR}, offset::OffsetCoord)
    q = offset.col - (offset.row - (offset.row & 1)) ÷ 2
    r = offset.row
    Axial(q, r)
end

function offset_to_axial(::Type{EvenR}, offset::OffsetCoord)
    q = offset.col - (offset.row + (offset.row & 1)) ÷ 2
    r = offset.row
    Axial(q, r)
end

function offset_to_axial(::Type{OddQ}, offset::OffsetCoord)
    q = offset.col
    r = offset.row - (offset.col - (offset.col & 1)) ÷ 2
    Axial(q, r)
end

function offset_to_axial(::Type{EvenQ}, offset::OffsetCoord)
    q = offset.col
    r = offset.row - (offset.col + (offset.col & 1)) ÷ 2
    Axial(q, r)
end

function axial_to_offset(::Type{OddR}, axial::Axial)
    col = axial.q + (axial.r - (axial.r & 1)) ÷ 2
    row = axial.r
    OffsetCoord(col, row)
end

function axial_to_offset(::Type{EvenR}, axial::Axial)
    col = axial.q + (axial.r + (axial.r & 1)) ÷ 2
    row = axial.r
    OffsetCoord(col, row)
end

function axial_to_offset(::Type{OddQ}, axial::Axial)
    col = axial.q
    row = axial.r + (axial.q - (axial.q & 1)) ÷ 2
    OffsetCoord(col, row)
end

function axial_to_offset(::Type{EvenQ}, axial::Axial)
    col = axial.q
    row = axial.r + (axial.q + (axial.q & 1)) ÷ 2
    OffsetCoord(col, row)
end

#-----------------------------------------------------------------------------# Arithmetic Operations
Base.:+(a::Cube, b::Cube) = Cube(a.q + b.q, a.r + b.r, a.s + b.s)
Base.:-(a::Cube, b::Cube) = Cube(a.q - b.q, a.r - b.r, a.s - b.s)
Base.:*(a::Cube, k::Real) = Cube(a.q * k, a.r * k, a.s * k)
Base.:*(k::Real, a::Cube) = a * k

Base.:+(a::Axial, b::Axial) = Axial(a.q + b.q, a.r + b.r)
Base.:-(a::Axial, b::Axial) = Axial(a.q - b.q, a.r - b.r)
Base.:*(a::Axial, k::Real) = Axial(a.q * k, a.r * k)
Base.:*(k::Real, a::Axial) = a * k

Base.:(==)(a::Cube, b::Cube) = a.q == b.q && a.r == b.r && a.s == b.s
Base.:(==)(a::Axial, b::Axial) = a.q == b.q && a.r == b.r
Base.:(==)(a::OffsetCoord, b::OffsetCoord) = a.col == b.col && a.row == b.row

Base.hash(a::Cube, h::UInt) = hash(a.q, hash(a.r, hash(a.s, h)))
Base.hash(a::Axial, h::UInt) = hash(a.q, hash(a.r, h))
Base.hash(a::OffsetCoord, h::UInt) = hash(a.col, hash(a.row, h))

#-----------------------------------------------------------------------------# Directions
# Direction vectors for flat-topped hexagons (pointy-top would be different)
const CUBE_DIRECTIONS = [
    Cube(1, 0, -1),   # East
    Cube(1, -1, 0),   # Northeast
    Cube(0, -1, 1),   # Northwest
    Cube(-1, 0, 1),   # West
    Cube(-1, 1, 0),   # Southwest
    Cube(0, 1, -1),   # Southeast
]

const AXIAL_DIRECTIONS = [
    Axial(1, 0),    # East
    Axial(1, -1),   # Northeast
    Axial(0, -1),   # Northwest
    Axial(-1, 0),   # West
    Axial(-1, 1),   # Southwest
    Axial(0, 1),    # Southeast
]

"""
    direction(coord, dir::Integer)

Get the direction vector for a given direction (1-6).
"""
direction(::Cube, dir::Integer) = CUBE_DIRECTIONS[dir]
direction(::Axial, dir::Integer) = AXIAL_DIRECTIONS[dir]

"""
    neighbor(coord, dir::Integer)

Get the neighboring hexagon in direction `dir` (1-6).
"""
neighbor(coord::Union{Cube,Axial}, dir::Integer) = coord + direction(coord, dir)

"""
    neighbors(coord)

Get all 6 neighboring hexagons.
"""
neighbors(coord::Cube) = [coord + d for d in CUBE_DIRECTIONS]
neighbors(coord::Axial) = [coord + d for d in AXIAL_DIRECTIONS]

#-----------------------------------------------------------------------------# Distance
"""
    distance(a::Cube, b::Cube)

Calculate the Manhattan distance between two hexagons in cube coordinates.
"""
distance(a::Cube, b::Cube) = (abs(a.q - b.q) + abs(a.r - b.r) + abs(a.s - b.s)) ÷ 2

"""
    distance(a::Axial, b::Axial)

Calculate the Manhattan distance between two hexagons in axial coordinates.
"""
function distance(a::Axial, b::Axial)
    (abs(a.q - b.q) + abs(a.r - b.r) + abs(a.q + a.r - b.q - b.r)) ÷ 2
end

#-----------------------------------------------------------------------------# Rounding
"""
    cube_round(cube::Cube)

Round floating-point cube coordinates to the nearest integer cube coordinates.
"""
function cube_round(cube::Cube)
    q = round(Int, cube.q)
    r = round(Int, cube.r)
    s = round(Int, cube.s)

    q_diff = abs(q - cube.q)
    r_diff = abs(r - cube.r)
    s_diff = abs(s - cube.s)

    if q_diff > r_diff && q_diff > s_diff
        q = -r - s
    elseif r_diff > s_diff
        r = -q - s
    else
        s = -q - r
    end

    Cube(q, r, s)
end

"""
    axial_round(axial::Axial)

Round floating-point axial coordinates to the nearest integer axial coordinates.
"""
axial_round(axial::Axial) = convert(Axial, cube_round(convert(Cube, axial)))

#-----------------------------------------------------------------------------# Linear Interpolation
"""
    lerp(a, b, t)

Linear interpolation between two values.
"""
lerp(a::Real, b::Real, t::Real) = a + (b - a) * t

"""
    cube_lerp(a::Cube, b::Cube, t::Real)

Linear interpolation between two cube coordinates.
"""
function cube_lerp(a::Cube, b::Cube, t::Real)
    Cube(lerp(a.q, b.q, t), lerp(a.r, b.r, t), lerp(a.s, b.s, t))
end

"""
    axial_lerp(a::Axial, b::Axial, t::Real)

Linear interpolation between two axial coordinates.
"""
function axial_lerp(a::Axial, b::Axial, t::Real)
    Axial(lerp(a.q, b.q, t), lerp(a.r, b.r, t))
end

#-----------------------------------------------------------------------------# Line Drawing
"""
    line(a::Cube, b::Cube)

Draw a line between two hexagons in cube coordinates. Returns a vector of hexagons.
"""
function line(a::Cube, b::Cube)
    N = distance(a, b)
    N == 0 && return [a]

    results = Vector{Cube{Int}}(undef, N + 1)
    for i in 0:N
        results[i + 1] = cube_round(cube_lerp(a, b, i / N))
    end
    results
end

"""
    line(a::Axial, b::Axial)

Draw a line between two hexagons in axial coordinates. Returns a vector of hexagons.
"""
function line(a::Axial, b::Axial)
    a_cube = convert(Cube, a)
    b_cube = convert(Cube, b)
    cube_line = line(a_cube, b_cube)
    [convert(Axial, c) for c in cube_line]
end

#-----------------------------------------------------------------------------# Range / Reachability
"""
    range(center, radius)

Get all hexagons within `radius` steps of `center`.
"""
function Base.range(center::Cube, radius::Integer)
    results = Cube{Int}[]
    for q in -radius:radius
        r1 = max(-radius, -q - radius)
        r2 = min(radius, -q + radius)
        for r in r1:r2
            push!(results, center + Cube(q, r, -q-r))
        end
    end
    results
end

function Base.range(center::Axial, radius::Integer)
    cube_results = range(convert(Cube, center), radius)
    [convert(Axial, c) for c in cube_results]
end

"""
    ring(center, radius)

Get all hexagons exactly `radius` steps from `center`.
"""
function ring(center::Cube, radius::Integer)
    radius == 0 && return [center]

    results = Cube{Int}[]
    hex = center + direction(center, 5) * radius

    for i in 1:6
        for j in 1:radius
            push!(results, hex)
            hex = neighbor(hex, i)
        end
    end
    results
end

function ring(center::Axial, radius::Integer)
    cube_results = ring(convert(Cube, center), radius)
    [convert(Axial, c) for c in cube_results]
end

"""
    spiral(center, radius)

Get all hexagons in a spiral pattern from `center` out to `radius`.
"""
function spiral(center::Union{Cube,Axial}, radius::Integer)
    results = [center]
    for k in 1:radius
        append!(results, ring(center, k))
    end
    results
end

#-----------------------------------------------------------------------------# Rotation
"""
    rotate_left(hex::Cube)

Rotate cube coordinates 60 degrees counter-clockwise around the origin.
"""
rotate_left(hex::Cube) = Cube(-hex.s, -hex.q, -hex.r)

"""
    rotate_right(hex::Cube)

Rotate cube coordinates 60 degrees clockwise around the origin.
"""
rotate_right(hex::Cube) = Cube(-hex.r, -hex.s, -hex.q)

"""
    rotate_left(hex::Axial)

Rotate axial coordinates 60 degrees counter-clockwise around the origin.
"""
function rotate_left(hex::Axial)
    cube = convert(Cube, hex)
    convert(Axial, rotate_left(cube))
end

"""
    rotate_right(hex::Axial)

Rotate axial coordinates 60 degrees clockwise around the origin.
"""
function rotate_right(hex::Axial)
    cube = convert(Cube, hex)
    convert(Axial, rotate_right(cube))
end

#-----------------------------------------------------------------------------# Reflection
"""
    reflect_q(hex::Cube)

Reflect across the Q axis.
"""
reflect_q(hex::Cube) = Cube(hex.q, hex.s, hex.r)

"""
    reflect_r(hex::Cube)

Reflect across the R axis.
"""
reflect_r(hex::Cube) = Cube(hex.s, hex.r, hex.q)

"""
    reflect_s(hex::Cube)

Reflect across the S axis.
"""
reflect_s(hex::Cube) = Cube(hex.r, hex.q, hex.s)

end # module Hexagons

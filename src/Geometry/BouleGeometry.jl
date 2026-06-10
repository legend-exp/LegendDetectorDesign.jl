# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

"""
    BouleGeometry{N, T} <: AbstractBouleGeometry{N, T}

Axially-symmetric crystalline boule profile: radius as a function of position
along the growth axis.

The boule is represented by `N` sample points `(zᵢ, radiusᵢ)` plus a
Fritsch–Carlson monotonic cubic Hermite interpolation through them, allowing
the radius to be queried at any axial coordinate via [`get_boule_radius`](@ref).
The interpolation supports non-uniform knots and is monotonicity-preserving —
it cannot overshoot above the local maximum or below the local minimum of the
samples, which matches the physics of a tapered boule. All lengths are stored
in internal length units (mm).

# Fields
- `z::SVector{N, T}`: axial sample positions, sorted ascending.
- `radius::SVector{N, T}`: boule radius at each `z`.
- `spline::Interpolations.MonotonicInterpolation{T}`: Fritsch–Carlson interpolation
  `radius(z)` used for evaluation.
"""
struct BouleGeometry{N, T} <: AbstractBouleGeometry{N, T}
    z::SVector{N, T}
    radius::SVector{N, T}
    spline::Interpolations.MonotonicInterpolation{T}
end

"""
    BouleGeometry(::Type{T}; z, radius) -> BouleGeometry{N, T}

Construct a [`BouleGeometry`](@ref) from measured `(z, radius)` samples.

`z` and `radius` may be `Number`s or `Vector`s of equal length, with or without
units (unitful inputs are converted to internal length units). If a single
scalar is given for `z`, the boule is treated as a uniform cylinder of length
`z` and constant `radius`. The samples are sorted by `z` and fit with
`Interpolations.FritschCarlsonMonotonicInterpolation` — a non-uniform-knot,
monotonicity-preserving cubic Hermite spline.

See also [`get_boule_radius`](@ref), [`get_physical_volume`](@ref).
"""
function BouleGeometry(::Type{T};
        z::Union{Number, Vector{<:Number}},
        radius::Union{Number, Vector{<:Number}}
    ) where {T}

    @assert length(z) == length(radius) "Vectors must be the same length"

    if length(z) == 1 z, radius = [0*unit(z), z], [radius, radius] end

    N = length(z)

    idx = sortperm(z)
    z      = SVector{N, T}(to_internal_length_units.(z[idx]))
    radius = SVector{N, T}(to_internal_length_units.(radius[idx]))

    spline = interpolate(z, radius, FritschCarlsonMonotonicInterpolation())
    BouleGeometry{N, T}(z, radius, spline)
end

"""
    get_boule_radius(geo::BouleGeometry, z::Number) -> Real

Interpolated boule radius at axial position `z`. `z` may be unitful or a bare
number in internal length units (mm); out-of-range queries are clamped to the
boule's `[z[1], z[end]]` (flat extrapolation). Returns the value in internal
length units.
"""
function get_boule_radius(geo::BouleGeometry, z::Number)
    z_internal = to_internal_length_units(z)
    geo.spline(clamp(z_internal, geo.z[1], geo.z[end]))
end

"""
    get_physical_volume(geo::BouleGeometry) -> Real

Return the physical volume of the boule crystal in cm³, computed analytically
as a sum of truncated-cone (frustum) volumes between adjacent `(z, radius)`
samples:

```math
V = \\sum_{i=1}^{N-1} \\frac{π}{3} (z_{i+1} - z_i)(r_i^2 + r_i r_{i+1} + r_{i+1}^2).
```

This is exact for a piecewise-linear radius profile, which is the
best interpretation of a finite set of measured radii.
"""
function get_physical_volume(geo::BouleGeometry{N,T})::T where {N,T}
    z, r = geo.z, geo.radius
    Vtot = (π/3) * sum((z[i+1] - z[i]) * (r[i]^2 + r[i]*r[i+1] + r[i+1]^2) for i in 1:N-1)
    ustrip(u"cm^3", Vtot * internal_length_unit^3)
end

function get_unicode_rep(::BouleGeometry; cut::Bool = false)
    if cut 
        " ╭──╮╭──────────╮╭───╮  ", 
        "/   ││          ││    \\ ",
        "│   ││          ││     )",
        "\\   ││          ││    / ",
        " ╰──╯╰──────────╯╰───╯  "
    else
        " ╭───────────────────╮  ", 
        "/                     \\ ",
        "│                      )",
        "\\                     / ",
        " ╰───────────────────╯  "
    end
end

function print(io::IO, geo::BouleGeometry)
    g1,g2,g3,g4,g5 = get_unicode_rep(geo)
    r, l = Int.(round.((maximum(geo.radius), geo.z[end]-geo.z[1])))
    println(io, "$g1  ╮    $(typeof(geo))")
    println(io, "$g2$(lpad(string(r), 3, ' '))ₘₘ  ╰─Radius: Maximum, Minimum")
    println(io, "$g3  ╯      ╰─$(maximum(geo.radius)), $(minimum(geo.radius)) mm")
    println(io, "$g4       ╰─Length: $(geo.z[end]-geo.z[1]) mm")
    println(io, "$g5       ╰─Mass: $(Int(round(get_physical_volume(geo)*ge_76_density))) g ")
    println(io, "╰────────$(lpad(string(l), 2, ' '))ₘₘ─────────╯ ")
end
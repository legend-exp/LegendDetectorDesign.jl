# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

struct BouleGeometry{N, T} <: AbstractBouleGeometry{N, T}
    z::SVector{N, T}
    radius::SVector{N, T}
    spline::Interpolations.MonotonicInterpolation{T}
end

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

function get_boule_radius(geo::BouleGeometry, z::Number)
    z_internal = to_internal_length_units(z)
    geo.spline(clamp(z_internal, geo.z[1], geo.z[end]))
end

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
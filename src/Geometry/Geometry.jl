# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

"""
    ValidGeometry

Type-level tag indicating that a geometry satisfies its internal consistency
checks (radii, heights, taper sizes). Used as the second type parameter of
[`AbstractDesignGeometry`](@ref) to enable dispatch on validity at compile time.

See also [`InvalidGeometry`](@ref) and [`is_valid_geometry`](@ref).
"""
abstract type ValidGeometry end

"""
    InvalidGeometry

Type-level tag indicating that a geometry failed at least one of its internal
consistency checks. Methods that require a physically meaningful geometry
(e.g. [`get_physical_volume`](@ref)) typically return `missing` or skip work
when dispatched on this tag.

See also [`ValidGeometry`](@ref) and [`is_valid_geometry`](@ref).
"""
abstract type InvalidGeometry end

"""
    AbstractDesignGeometry{T <: SSDFloat, V <: Union{ValidGeometry, InvalidGeometry}}

Supertype for detector design geometries (e.g. [`InvertedCoaxGeometry`](@ref)).

`T` is the floating-point precision used for dimensions, and `V` is a validity
tag — either [`ValidGeometry`](@ref) or [`InvalidGeometry`](@ref) — set by the
constructor from simple consistency checks on the inputs.
"""
abstract type AbstractDesignGeometry{T <: SSDFloat, V <: Union{ValidGeometry, InvalidGeometry}} end

"""
    AbstractBouleGeometry{N, T <: SSDFloat}

Supertype for crystalline-boule geometries (e.g. [`BouleGeometry`](@ref)).
Parameters: `N` is the number of axial sample points, `T` the float type.
"""
abstract type AbstractBouleGeometry{N, T <: SSDFloat} end

include("InvertedCoaxGeometry.jl")
include("BouleGeometry.jl")

"""
    get_physical_volume(geo::AbstractDesignGeometry{T, InvalidGeometry}) -> Missing

Fallback returning `missing` for geometries tagged [`InvalidGeometry`](@ref).
The volume of a non-physical geometry is not well defined, so downstream
calculations that consume the result (e.g. mass) propagate `missing` cleanly
instead of producing a misleading number.
"""
get_physical_volume(geo::AbstractDesignGeometry{T, InvalidGeometry}) where {T} = missing

"""
    is_valid_geometry(geo::AbstractDesignGeometry) -> Bool

`true` if `geo` is tagged [`ValidGeometry`](@ref), `false` if [`InvalidGeometry`](@ref).
Equivalent to a check on the second type parameter; provided for readability at
call sites.
"""
is_valid_geometry(geo::AbstractDesignGeometry{T, InvalidGeometry}) where {T} = false
is_valid_geometry(geo::AbstractDesignGeometry{T, ValidGeometry}) where {T} = true

"""
    meta_to_geo(::Type{T}, meta::PropDict) -> AbstractDesignGeometry{T}

Dispatch on `meta.type` and construct the matching design geometry with element
type `T` from a LEGEND detector-metadata `PropDict`. Currently supports
`"icpc"` (→ [`InvertedCoaxGeometry`](@ref)); throws `ArgumentError` for any
other `meta.type`.
"""
function meta_to_geo(::Type{T}, meta::PropDict) where {T}
    if meta.type == "icpc"
        meta_to_geo(InvertedCoaxGeometry{T}, meta)
    else
        throw(ArgumentError("Unknown geometry type $(meta.type)"))
    end
end

function show(io::IO, geo::Union{AbstractBouleGeometry, AbstractDesignGeometry{<:Any, ValidGeometry}})
    print(io, geo)
end

function show(io::IO, ::MIME"text/plain", geo::Union{AbstractBouleGeometry, AbstractDesignGeometry{<:Any, ValidGeometry}})
    show(io, geo)
end
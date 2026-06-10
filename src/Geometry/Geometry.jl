# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

abstract type ValidGeometry end
abstract type InvalidGeometry end
abstract type AbstractDesignGeometry{T <: SSDFloat, V <: Union{ValidGeometry, InvalidGeometry}} end
abstract type AbstractBouleGeometry{N, T <: SSDFloat} end

include("InvertedCoaxGeometry.jl")
include("BouleGeometry.jl")

get_physical_volume(geo::AbstractDesignGeometry{T, InvalidGeometry}) where {T} = missing

is_valid_geometry(geo::AbstractDesignGeometry{T, InvalidGeometry}) where {T} = false
is_valid_geometry(geo::AbstractDesignGeometry{T, ValidGeometry}) where {T} = true

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
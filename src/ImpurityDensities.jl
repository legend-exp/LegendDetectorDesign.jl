# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

struct SegregationImpurityDensity{T} <: AbstractImpurityDensity{T}
    ρ0::T
    ρm::T
    ρl::T
    l::T
    ml::T
end

function SolidStateDetectors.get_impurity_density(segdm::SegregationImpurityDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T}
    pt::CartesianPoint{T} = CartesianPoint(pt)
    return T(segdm.ρ0 + segdm.ρm*z + segdm.ρl*exp((z-segdm.l)/segdm.ml))
end
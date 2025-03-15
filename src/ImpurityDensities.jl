# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

@. segcurve(x,p) = p[1] + p[2]*x + p[3]*exp((x-p[4])/p[5])

struct SegregationImpurityDensity{T} <: AbstractImpurityDensity{T}
    # a + b*z + c*exp((z-L)/tau) -> needs at least 4 points
    a::T
    b::T 
    c::T 
    L::T
    tau::T
    det_z0::T
end

(*)(idm::SegregationImpurityDensity{T}, scale::Real) where {T} = (*)(scale, idm)

(*)(scale::Real, idm::SegregationImpurityDensity{T}) where {T} = SegregationImpurityDensity{T}(T(scale*idm.a), T(scale*idm.b), T(scale*idm.c), idm.L, idm.tau, idm.det_z0)

function SegregationImpurityDensity{T}(segfitpars::Vector{<:Quantity}, det_z0::Number) where {T}
    @assert length(segfitpars) == 5 "Need 5 fit parameters"
    SegregationImpurityDensity(
        T.(to_internal_ssd_units.(segfitpars))..., 
        T(to_internal_ssd_units(internal_length_unit*to_internal_length_units(det_z0)))
        )
end

function SegregationImpurityDensity{T}(segfitpars::Vector{T}, det_z0::T) where {T}
    @assert length(segfitpars) == 5 "Need 5 fit parameters"
    SegregationImpurityDensity(
        T.(to_internal_ssd_units.(internal_segregation_fit_units .* segfitpars))...,
        T(to_internal_ssd_units(internal_length_unit*det_z0))
        )
end

function SolidStateDetectors.get_impurity_density(
    idm::SegregationImpurityDensity, pt::AbstractCoordinatePoint{T}
    )::T where {T}
    cpt = CartesianPoint(pt)
    z = cpt[3]

    # the function parameters are in crystal axis coordinates i.e. z = 0 is seed end, z = L crystal length 
    # -> convert to detector coordiantes where z = 0 corresponds to p+ contact i.e. z -> det_z0 - z
    idm.a .+ idm.b * (idm.det_z0 .- z) .+ idm.c * exp.((idm.det_z0 .- z .- idm.L)/idm.tau)
end
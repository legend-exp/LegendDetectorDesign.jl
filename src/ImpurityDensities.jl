# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

@. segcurve(z,p) = p[1] + p[2]*z + p[3]*exp((z-p[4])/p[5])
@. lin(z,p) = p[1] + p[2]*z

#=struct SegregationImpurityDensity{T} <: AbstractImpurityDensity{T}
    # a + b*z + c*exp((z-L)/tau) -> needs at least 4 points
    a::T
    b::T 
    c::T 
    L::T
    tau::T
    det_z0::T
end
=#
const internal_segregation_fit_units = [internal_impurity_quantity, 
    internal_impurity_quantity/internal_length_unit, 
    internal_impurity_quantity, 
    1*internal_length_unit, 
    1*internal_length_unit]
#=
struct LinearImpurityDensity{T} <: AbstractImpurityDensity{T}
    # a + b*z
    a::T
    b::T 
    det_z0::T
end
=#
const internal_linear_fit_units = [internal_impurity_quantity, internal_impurity_quantity/internal_length_unit]

(*)(idm::AbstractImpurityDensity{T}, scale::Real) where {T} = (*)(scale, idm)

(*)(scale::Real, idm::LinBouleImpurityDensity{T}) where {T} = LinBouleImpurityDensity{T}(T(scale*idm.a), T(scale*idm.b), idm.det_z0)

(*)(scale::Real, idm::ParBouleImpurityDensity{T}) where {T} = ParBouleImpurityDensity{T}(T(scale*idm.a), T(scale*idm.b), T(scale*idm.c), idm.det_z0)

(*)(scale::Real, idm::LinExpBouleImpurityDensity{T}) where {T} = LinExpBouleImpurityDensity{T}(T(scale*idm.a), T(scale*idm.b), T(scale*idm.n), idm.l, idm.m, idm.det_z0)

(*)(scale::Real, idm::ParExpBouleImpurityDensity{T}) where {T} = ParExpBouleImpurityDensity{T}(T(scale*idm.a), T(scale*idm.b), T(scale*idm.c), T(scale*idm.n), idm.l, idm.m, idm.det_z0)


function LinBouleImpurityDensity{T}(fitpars::Vector{<:Quantity}, det_z0::Number) where {T}
    @assert length(fitpars) == 2 "Need 2 fit parameters"
    LinBouleImpurityDensity(
        T.(to_internal_ssd_units.(fitpars))..., 
        T(to_internal_ssd_units(internal_length_unit*to_internal_length_units(det_z0)))
        )
end

function LinBouleImpurityDensity{T}(fitpars::Vector{T}, det_z0::T) where {T}
    @assert length(fitpars) == 2 "Need 2 fit parameters"
    LinBouleImpurityDensity(
        T.(to_internal_ssd_units.(internal_linear_fit_units .* fitpars))...,
        T(to_internal_ssd_units(internal_length_unit*det_z0))
        )
end

function ParBouleImpurityDensity{T}(fitpars::Vector{<:Quantity}, det_z0::Number) where {T}
    @assert length(fitpars) == 3 "Need 3 fit parameters"
    ParBouleImpurityDensity(
        T.(to_internal_ssd_units.(fitpars))..., 
        T(to_internal_ssd_units(internal_length_unit*to_internal_length_units(det_z0)))
        )
end

function ParBouleImpurityDensity{T}(fitpars::Vector{T}, det_z0::T) where {T}
    @assert length(fitpars) == 3 "Need 3 fit parameters"
    ParBouleImpurityDensity(
        T.(to_internal_ssd_units.(internal_linear_fit_units .* fitpars))...,
        T(to_internal_ssd_units(internal_length_unit*det_z0))
        )
end

function LinExpBouleImpurityDensity{T}(fitpars::Vector{<:Quantity}, det_z0::Number) where {T}
    @assert length(fitpars) == 5 "Need 5 fit parameters"
    LinExpBouleImpurityDensity(
        T.(to_internal_ssd_units.(fitpars))..., 
        T(to_internal_ssd_units(internal_length_unit*to_internal_length_units(det_z0)))
        )
end

function LinExpBouleImpurityDensity{T}(fitpars::Vector{T}, det_z0::T) where {T}
    @assert length(fitpars) == 5 "Need 5 fit parameters"
    LinExpBouleImpurityDensity(
        T.(to_internal_ssd_units.(internal_segregation_fit_units .* fitpars))...,
        T(to_internal_ssd_units(internal_length_unit*det_z0))
        )
end

function ParExpBouleImpurityDensity{T}(fitpars::Vector{<:Quantity}, det_z0::Number) where {T}
    @assert length(fitpars) == 5 "Need 6 fit parameters"
    ParExpBouleImpurityDensity(
        T.(to_internal_ssd_units.(fitpars))..., 
        T(to_internal_ssd_units(internal_length_unit*to_internal_length_units(det_z0)))
        )
end

function ParExpBouleImpurityDensity{T}(fitpars::Vector{T}, det_z0::T) where {T}
    @assert length(fitpars) == 5 "Need 6 fit parameters"
    ParExpBouleImpurityDensity(
        T.(to_internal_ssd_units.(internal_segregation_fit_units .* fitpars))...,
        T(to_internal_ssd_units(internal_length_unit*det_z0))
        )
end
#=
function SolidStateDetectors.get_impurity_density(
    idm::SegregationImpurityDensity, pt::AbstractCoordinatePoint{T}
    )::T where {T}
    cpt = CartesianPoint(pt)
    z = cpt[3]

    # the function parameters are in crystal axis coordinates i.e. z = 0 is seed end, z = L crystal length 
    # -> convert to detector coordiantes where z = 0 corresponds to p+ contact i.e. z -> det_z0 - z
    idm.a .+ idm.b * (idm.det_z0 .- z) .+ idm.c * exp.((idm.det_z0 .- z .- idm.L)/idm.tau)
end

function SolidStateDetectors.get_impurity_density(
    idm::LinearImpurityDensity, pt::AbstractCoordinatePoint{T}
    )::T where {T}
    cpt = CartesianPoint(pt)
    z = cpt[3]

    # the function parameters are in crystal axis coordinates i.e. z = 0 is seed end, z = L crystal length 
    # -> convert to detector coordiantes where z = 0 corresponds to p+ contact i.e. z -> det_z0 - z
    idm.a .+ idm.b * (idm.det_z0 .- z)
end
=#
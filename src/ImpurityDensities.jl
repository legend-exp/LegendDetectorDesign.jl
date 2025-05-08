# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

function get_default_xtal_meta(det::DetectorDesign)
    PropDict(
        :impurity_curve => PropDict(
                :model => "constant_boule",
                :parameters => PropDict("value" => 0,)
                ),
        :slices => PropDict(Symbol(det.name[end]) => PropDict("detector_offset_in_mm" => 0))
    )
end

fit_function(model::Symbol) = fit_function(Val(model))
from_internal_units(model::Symbol, impurity_unit::Unitful.Units, length_unit::Unitful.Units, x::Vector{<:Real}) =  from_internal_units(Val(model), impurity_unit, length_unit, x)

@. linear_boule(z,p) = p[1] + p[2]*z
fit_function(::Val{:linear_boule}) = linear_boule

function from_internal_units(::Val{:linear_boule}, impurity_unit::Unitful.Units, length_unit::Unitful.Units, x::Vector{<:Real})
    [
        uconvert(impurity_unit, x[1] * internal_impurity_quantity),
        uconvert(impurity_unit / length_unit, x[2] * internal_impurity_quantity / internal_length_unit),
    ]
end

@. parabolic_boule(z,p) = p[1] + p[2]*z + p[3]*z^2
fit_function(::Val{:parabolic_boule}) = parabolic_boule

function from_internal_units(::Val{:parabolic_boule}, impurity_unit::Unitful.Units, length_unit::Unitful.Units, x::Vector{<:Real})
    [
        uconvert(impurity_unit, x[1] * internal_impurity_quantity),
        uconvert(impurity_unit / length_unit, x[2] * internal_impurity_quantity / internal_length_unit),
        uconvert(impurity_unit / length_unit^2, x[3] * internal_impurity_quantity / internal_length_unit^2),
    ]
end

@. linear_exponential_boule(z,p) = p[1] + p[2]*z + p[3]*exp((z-p[4])/p[5])
fit_function(::Val{:linear_exponential_boule}) = linear_exponential_boule

function from_internal_units(::Val{:linear_exponential_boule}, impurity_unit::Unitful.Units, length_unit::Unitful.Units, x::Vector{<:Real})
    [
        uconvert(impurity_unit, x[1] * internal_impurity_quantity),
        uconvert(impurity_unit / length_unit, x[2] * internal_impurity_quantity / internal_length_unit),
        uconvert(impurity_unit, x[3] * internal_impurity_quantity),
        uconvert(length_unit, x[4] * internal_length_unit),
        uconvert(length_unit, x[5] * internal_length_unit)
    ]
end
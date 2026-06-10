# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

"""
    get_default_xtal_meta(det::DetectorDesign) -> PropDict

Build a minimal crystal-metadata `PropDict` for `det`, suitable as the
`xtal_meta` argument to SolidStateDetectors' `Simulation` and
`SolidStateDetector` constructors.

The impurity model is set to `"constant_boule"` with value `0` — a placeholder
that SSD parses but does not use, because LegendDetectorDesign supplies the
real impurity density via a separate `AbstractImpurityDensity` object passed to
the SSD constructor.

The `:slices` entry encodes the boule order letter (the last character of
`det.name`) with a zero detector offset.
"""
function get_default_xtal_meta(det::DetectorDesign)
    PropDict(
        :impurity_curve => PropDict(
                :model => "constant_boule",
                :parameters => PropDict("value" => 0,)
                ),
        :slices => PropDict(Symbol(det.name[end]) => PropDict("detector_offset_in_mm" => 0))
    )
end

"""
    fit_function(model::Symbol) -> Function
    impurity_density_model(model::Symbol) -> Type{<:AbstractImpurityDensity}
    fit_parameter_names(model::Symbol) -> Vector{Symbol}
    fit_parameter_units(model::Symbol) -> Vector{Unitful.Units}

Look up the four pieces of metadata for a named boule impurity model by symbol:

- `fit_function` — the numerical curve `ρ(z, p)` (broadcasted over `z`) used to
  fit measured impurity data and to evaluate the impurity density at a point.
- `impurity_density_model` — the SSD `AbstractImpurityDensity` subtype that
  realises this model inside a simulation.
- `fit_parameter_names` / `fit_parameter_units` — parameter names and the
  units expected for each parameter.

Supported `model` symbols: `:linear_boule`, `:parabolic_boule`,
`:linear_exponential_boule`, `:parabolic_exponential_boule`. Each dispatches
through `Val(model)`.
"""
fit_function(model::Symbol) = fit_function(Val(model))
impurity_density_model(model::Symbol) = impurity_density_model(Val(model))
fit_parameter_names(model::Symbol) =  fit_parameter_names(Val(model))
fit_parameter_units(model::Symbol) =  fit_parameter_units(Val(model))

"""
    linear_boule(z, p) = p[1] + p[2]*z

Two-parameter linear impurity profile `ρ(z) = a + b·z` (broadcasted over `z`).
Parameter order: `a, b` in [`internal_impurity_quantity`](@ref) and
`impurity / mm` respectively.
"""
@. linear_boule(z,p) = p[1] + p[2]*z
fit_function(::Val{:linear_boule}) = linear_boule
impurity_density_model(::Val{:linear_boule}) = LinBouleImpurityDensity
fit_parameter_names(::Val{:linear_boule}) = Symbol.(["a", "b"])

function fit_parameter_units(::Val{:linear_boule})
    [
        internal_impurity_quantity,
        internal_impurity_quantity / internal_length_unit
    ]
end

"""
    parabolic_boule(z, p) = p[1] + p[2]*z + p[3]*z^2

Three-parameter quadratic impurity profile (`a, b, c`).
"""
@. parabolic_boule(z,p) = p[1] + p[2]*z + p[3]*z^2
fit_function(::Val{:parabolic_boule}) = parabolic_boule
impurity_density_model(::Val{:parabolic_boule}) = ParBouleImpurityDensity
fit_parameter_names(::Val{:parabolic_boule}) = Symbol.(["a", "b", "c"])

function fit_parameter_units(::Val{:parabolic_boule})
    [
        internal_impurity_quantity,
        internal_impurity_quantity / internal_length_unit,
        internal_impurity_quantity / internal_length_unit^2
    ]
end

"""
    linear_exponential_boule(z, p) = p[1] + p[2]*z + p[3]*exp((z - p[4]) / p[5])

Five-parameter "linear + exponential tail" impurity profile. Parameters
`(a, b, n, l, m)`: linear baseline `a + b·z` plus an exponential rise with
amplitude `n`, onset position `l`, and decay length `m`. Typical fit form for
boules with a sharp tail-end gradient at the crystal cut.
"""
@. linear_exponential_boule(z,p) = p[1] + p[2]*z + p[3]*exp((z-p[4])/p[5])
fit_function(::Val{:linear_exponential_boule}) = linear_exponential_boule
impurity_density_model(::Val{:linear_exponential_boule}) = LinExpBouleImpurityDensity
fit_parameter_names(::Val{:linear_exponential_boule}) = Symbol.(["a", "b", "n", "l", "m"])

function fit_parameter_units(::Val{:linear_exponential_boule})
    [
        internal_impurity_quantity,
        internal_impurity_quantity / internal_length_unit,
        internal_impurity_quantity,
        internal_length_unit,
        internal_length_unit
    ]
end

"""
    parabolic_exponential_boule(z, p) = p[1] + p[2]*z + p[3]*z^2 + p[4]*exp((z - p[5]) / p[6])

Six-parameter "quadratic + exponential tail" impurity profile. Parameters
`(a, b, c, n, l, m)`: quadratic baseline plus the same exponential tail as
[`linear_exponential_boule`](@ref). Use when the bulk profile has noticeable
curvature in addition to the tail.
"""
@. parabolic_exponential_boule(z,p) = p[1] + p[2]*z + p[3]*z^2 + p[4]*exp((z-p[5])/p[6])
fit_function(::Val{:parabolic_exponential_boule}) = parabolic_exponential_boule
impurity_density_model(::Val{:parabolic_exponential_boule}) = ParExpBouleImpurityDensity
fit_parameter_names(::Val{:parabolic_exponential_boule}) = Symbol.(["a", "b", "c", "n", "l", "m"])

function fit_parameter_units(::Val{:parabolic_exponential_boule})
    [
        internal_impurity_quantity,
        internal_impurity_quantity / internal_length_unit,
        internal_impurity_quantity / internal_length_unit^2,
        internal_impurity_quantity,
        internal_length_unit,
        internal_length_unit
    ]
end
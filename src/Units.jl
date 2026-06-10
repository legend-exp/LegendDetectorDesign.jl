# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

"""
    ssd_ptype = -1

Sign convention used when converting impurity densities to
SolidStateDetectors internal units. SSD expresses impurity densities with a
sign that flips for p-type vs n-type material; LEGEND detectors are p-type by
default, so this multiplier is `-1`.
"""
const ssd_ptype = -1

"""
Module-level unit conventions for LegendDetectorDesign. All bare numbers stored
in geometry / design / boule structs are interpreted in these units:

| constant                    | value         | quantity                |
|-----------------------------|---------------|-------------------------|
| `internal_length_unit`      | `u"mm"`       | length                  |
| `internal_mass_unit`        | `u"g"`        | mass                    |
| `internal_angle_unit`       | `u"°"`        | angle                   |
| `internal_voltage_unit`     | `u"V"`        | voltage                 |
| `internal_efield_unit`      | `u"V/cm"`     | electric field          |
| `internal_density_unit`     | `u"g/cm^3"`   | mass density            |
| `internal_impurity_quantity`| `1e9u"cm^-3"` | impurity number density |

Use [`to_internal_units`](@ref) / [`from_internal_units`](@ref) to round-trip
through these conventions.
"""
const internal_length_unit  = u"mm"
const internal_mass_unit = u"g"
const internal_angle_unit   = u"°"
const internal_voltage_unit = u"V"
const internal_efield_unit  = u"V/cm"
const internal_density_unit  = u"g/cm^3"
const internal_impurity_quantity  = 1e9u"cm^-3"

"""
    ge_76_density = 5.544

Mass density of enriched ⁷⁶Ge in g/cm³, used to convert geometric volumes into
detector / boule mass.
"""
const ge_76_density = 5.544

"""
    default_operational_V = 5000

Fallback recommended operating voltage in volts, used when no `Vop` is supplied
when serialising a design to metadata.
"""
const default_operational_V = 5000

"""
    to_internal_units(x::Quantity) -> Real
    to_internal_units(::Missing) -> Missing

Strip the units from `x` and return its magnitude in the package's internal
units (see [`internal_length_unit`](@ref) and friends). Dispatches on the
dimension of `x`; throws `ArgumentError` for any dimension this package does
not handle (voltage, electric field, length, angle, mass, density, impurity
density). `missing` passes through unchanged so unitful broadcasts over
optionally-missing vectors keep working.
"""
to_internal_units(x::Quantity{<:Real}) = throw(ArgumentError("Unit $(unit(x)) unknown to LegendGeDesign.jl"))
to_internal_units(x::Quantity{<:Real, dimension(internal_voltage_unit)}) = ustrip(internal_voltage_unit, x)
to_internal_units(x::Quantity{<:Real, dimension(internal_efield_unit)})  = ustrip(internal_efield_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_length_unit)})  = ustrip(internal_length_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_angle_unit)})  = ustrip(internal_angle_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_mass_unit)})  = ustrip(internal_mass_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_density_unit)})  = ustrip(internal_density_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_impurity_quantity)})  = ustrip(x/internal_impurity_quantity)
to_internal_units(::Missing) = missing

"""
    to_internal_length_units(x) -> Real

Convert `x` to a bare number in `internal_length_unit` (mm). Accepts a `Real`
(returned as-is) or a `Quantity` with length dimension. Convenience wrapper
that avoids the dimension-check dispatch overhead of [`to_internal_units`](@ref)
when the dimension is already known to be length.
"""
to_internal_length_units(x::Real) = x
to_internal_length_units(x::Quantity{<:Real, dimension(internal_length_unit)}) = to_internal_units(x)

"""
    to_internal_angle_units(x) -> Real

Convert `x` to a bare number in `internal_angle_unit` (°). Accepts a `Real`
or an angle `Quantity`. Length-only counterpart of [`to_internal_length_units`](@ref).
"""
to_internal_angle_units(x::Real) = x
to_internal_angle_units(x::Quantity{<:Real, dimension(internal_angle_unit)}) = to_internal_units(x)

"""
    to_internal_ssd_units(x::Quantity) -> Real

Convert `x` to a bare number in **SolidStateDetectors**' internal units
(meters and m⁻³), not this package's. Used when handing dimensions or impurity
densities/gradients to SSD constructors. Impurity-density quantities are also
multiplied by [`ssd_ptype`](@ref) to match SSD's sign convention.
"""
to_internal_ssd_units(x::Quantity{<:Real, dimension(internal_length_unit)}) = ustrip(SolidStateDetectors.internal_length_unit,  x)
to_internal_ssd_units(x::Quantity{<:Real, dimension(internal_impurity_quantity)}) = ssd_ptype * ustrip(SolidStateDetectors.internal_length_unit^-3,  x)
to_internal_ssd_units(x::Quantity{<:Real, dimension(internal_impurity_quantity/internal_length_unit)}) = ssd_ptype * ustrip(SolidStateDetectors.internal_length_unit^-4,  x)
to_internal_ssd_units(x::Quantity{<:Real, dimension(internal_impurity_quantity/internal_length_unit^2)}) = ssd_ptype * ustrip(SolidStateDetectors.internal_length_unit^-5,  x)


"""
    from_internal_units(x::Real, unit) -> Quantity

Re-attach units to a bare number `x` previously stripped via
[`to_internal_units`](@ref), and convert it to the requested `unit`. Currently
defined for length dimensions only.
"""
from_internal_units(x::Real, unit::Unitful.Units{<:Any, dimension(internal_length_unit)}) = uconvert(unit, x * internal_length_unit);
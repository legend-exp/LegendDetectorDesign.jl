# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

const internal_length_unit  = u"mm"
const internal_mass_unit = u"g"
const internal_angle_unit   = u"°"
const internal_voltage_unit = u"V"
const internal_efield_unit  = u"V/cm"
const internal_impurity_unit  = u"1e9/cm^3"
const internal_density_unit  = u"g/cm^3"
const ge_76_density = 5.544
const default_operational_V = 5000

to_internal_units(x::Quantity{<:Real}) = throw(ArgumentError("Unit $(unit(x)) unknown to LegendGeDesign.jl"))
to_internal_units(x::Quantity{<:Real, dimension(internal_voltage_unit)}) = ustrip(internal_voltage_unit, x)
to_internal_units(x::Quantity{<:Real, dimension(internal_efield_unit)})  = ustrip(internal_efield_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_length_unit)})  = ustrip(internal_length_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_angle_unit)})  = ustrip(internal_angle_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_mass_unit)})  = ustrip(internal_mass_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_density_unit)})  = ustrip(internal_density_unit,  x)

to_internal_length_units(x::Real) = x
to_internal_length_units(x::Quantity{<:Real, dimension(internal_length_unit)}) = to_internal_units(x)
to_internal_angle_units(x::Real) = x
to_internal_angle_units(x::Quantity{<:Real, dimension(internal_angle_unit)}) = to_internal_units(x)

from_internal_units(x::Real, unit::Unitful.Units{<:Any, dimension(internal_length_unit)}) = uconvert(unit, x * internal_length_unit);
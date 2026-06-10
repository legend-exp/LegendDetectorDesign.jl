# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

"""
    AbstractCrystallineBoule{T <: SSDFloat}

Supertype for grown crystal boules carrying impurity / resistivity / Hall data
along their axis (e.g. [`CrystallineBoule`](@ref)).
"""
abstract type AbstractCrystallineBoule{T <: SSDFloat} end

"""
    CrystallineBoule{T} <: AbstractCrystallineBoule{T}

Container for everything known about a single grown germanium boule: its
external profile, the impurity-density model parameters fit from
manufacturer measurements, and the raw resistivity / Hall measurements those
parameters were fit to.

Detectors are cut from a boule at specific axial positions; coupling a
`CrystallineBoule` with a [`DetectorDesign`](@ref) provides the impurity
density inside that detector's volume.

# Fields
- `name::AbstractString`: boule identifier.
- `order::AbstractString`: ordering letter / code identifying the boule batch.
- `impurity_model::Function`: numerical curve `ρ(z, p)` — one of
  [`linear_boule`](@ref), [`parabolic_boule`](@ref),
  [`linear_exponential_boule`](@ref), [`parabolic_exponential_boule`](@ref).
- `geometry::Union{BouleGeometry, Missing}`: external boule profile.
- `impurity_resistivity`, `z_resistivity`: raw resistivity samples and their
  axial positions (both in internal units).
- `impurity_hall`, `z_hall`: Hall-measured impurity densities and their axial
  positions.
- `impurity_model_parameters`: best-fit parameters of `impurity_model`, in the
  order returned by `fit_parameter_names`.
"""
mutable struct CrystallineBoule{T} <: AbstractCrystallineBoule{T}
    name::AbstractString
    order::AbstractString
    impurity_model::Function
    geometry::Union{BouleGeometry{<:Any,T}, Missing}
    impurity_resistivity::Union{Vector{T}, Missing}
    z_resistivity::Union{Vector{T}, Missing}
    impurity_hall::Union{Vector{T}, Missing}
    z_hall::Union{Vector{T}, Missing}
    impurity_model_parameters::Union{Vector{T}, Missing}
end

"""
    CrystallineBoule(::Type{T};
                     name, order, impurity_model,
                     geometry = missing,
                     impurity_resistivity = missing, z_resistivity = missing,
                     impurity_hall = missing, z_hall = missing,
                     impurity_model_parameters = missing) -> CrystallineBoule{T}

Construct a [`CrystallineBoule`](@ref) with element type `T`.

`impurity_model` is a `Symbol` naming a registered curve (see
[`fit_function`](@ref)); the corresponding `impurity_model_parameters` must be
ordered to match `fit_parameter_names` for that model. All vector
inputs may be unitful — they are converted to internal units on construction
and stored as bare numbers. Any of the measurement vectors may be `missing`
when the data is not available.
"""
function CrystallineBoule(::Type{T};
        name::AbstractString,
        order::AbstractString,
        impurity_model::Symbol,
        geometry::Union{BouleGeometry{<:Any, T}, Missing} = missing,
        impurity_resistivity::Union{Vector{<:Number}, Missing} = missing,
        z_resistivity::Union{Vector{<:Number}, Missing} = missing,
        impurity_hall::Union{Vector{<:Number}, Missing} = missing,
        z_hall::Union{Vector{<:Number}, Missing} = missing,
        impurity_model_parameters::Union{Vector{<:Number}, Missing} = missing
    ) where {T <: SSDFloat}
    CrystallineBoule{T}(name, order, fit_function(impurity_model), geometry, to_internal_units.(impurity_resistivity), to_internal_units.(z_resistivity), to_internal_units.(impurity_hall), to_internal_units.(z_hall), to_internal_units.(impurity_model_parameters))
end

"""
    get_unitful_property(boule::CrystallineBoule, prop::Symbol)

Return a field of `boule` re-attached with its physical units, via `Val`
dispatch on `prop`. Currently supports `:impurity_model_parameters`, which
attaches the per-parameter units returned by `fit_parameter_units` for
the boule's impurity model.
"""
get_unitful_property(boule::CrystallineBoule, prop::Symbol) = get_unitful_property(boule, Val(prop))

get_unitful_property(boule::CrystallineBoule, ::Val{:impurity_model_parameters}) = boule.impurity_model_parameters .* fit_parameter_units(nameof(boule.impurity_model))

"""
    get_impurity_density(boule::CrystallineBoule, z) -> Quantity

Evaluate the boule's impurity-density model at axial position `z` and return
the value with units of `internal_impurity_quantity` (`1e9·cm⁻³`).
`z` is converted to internal units; vector / array `z` returns an array via
broadcasting.
"""
get_impurity_density(boule::CrystallineBoule, z::Number) = boule.impurity_model(to_internal_units(z), boule.impurity_model_parameters) * internal_impurity_quantity

get_impurity_density(boule, z::AbstractArray) = broadcast(z -> get_impurity_density(boule, z), z)

function print(io::IO, boule::CrystallineBoule)
    geo = boule.geometry
    g1,g2,g3,g4,g5 = get_unicode_rep(geo, cut = !ismissing(boule.z_hall))
    r, l = Int.(round.((maximum(geo.radius), geo.z[end]-geo.z[1])))
    println(io, "$g1  $(typeof(boule)) - $(boule.name)")
    println(io, "$g2  ╰─Impurity model: $(boule.impurity_model)")
    println(io, "$g3    ╰─Params: $(boule.impurity_model_parameters)")
    println(io, "$g4  ╰─Length: $(geo.z[end]-geo.z[1]) mm")
    println(io, "$g5  ╰─Mass: $(Int(round(get_physical_volume(geo)*ge_76_density))) g ")
end

function show(io::IO, boule::CrystallineBoule)
    print(io, boule)
end

function show(io::IO, ::MIME"text/plain", boule::CrystallineBoule)
    show(io, boule)
end

"""
    boule_to_meta(boule::CrystallineBoule, det::DetectorDesign) -> OrderedDict

Serialise `boule` together with the cutting offset of `det` into the
LEGEND-format boule-metadata `OrderedDict`, matching the layout used by the
LEGEND metadata YAMLs.

The returned dict contains:
- the boule name (last three characters of `boule.name`) and `order`,
- the raw Hall measurements (`impurity_measurements`: value in `1e9 cm⁻³` and
  axial distance from the seed end in mm),
- the impurity-curve model + parameters rounded to 4 significant figures,
- a `slices` entry mapping the detector's order letter to its axial offset
  from the seed end.

Companion of [`geo_to_meta`](@ref) and [`design_to_meta`](@ref). To write the
result to disk, call `YAML.write_file(path, boule_to_meta(boule, det))`.
"""
function boule_to_meta(boule::CrystallineBoule, det::DetectorDesign)
    OrderedDict(
        :name => boule.name[end-2:end],
        :order => boule.order,
        :impurity_measurements => OrderedDict(
                :value_in_1e9e_cm3 => boule.impurity_hall,
                :distance_from_seed_end_mm => boule.z_hall
        ),
        :impurity_curve => OrderedDict(
                :model => Symbol(boule.impurity_model),
                :parameters => OrderedDict(zip(fit_parameter_names(nameof(boule.impurity_model)), round.(boule.impurity_model_parameters, sigdigits = 4)))
                ),
        :slices => OrderedDict(
            Symbol(det.name[end]) => OrderedDict(
                :detector_offset_in_mm => det.offset
            )
        )
    )
end

function SolidStateDetectors.Simulation{T}(det::DetectorDesign{T}, boule::CrystallineBoule{T}, env::HPGeEnvironment = HPGeEnvironment(); kwargs...) where {T<:AbstractFloat}
    imp_model = ssd_ptype*impurity_density_model(nameof(boule.impurity_model)){T}(get_unitful_property(boule, :impurity_model_parameters), get_unitful_property(det, :offset))
    Simulation{T}(det, imp_model, env; kwargs...)
end
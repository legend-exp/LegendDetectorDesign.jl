# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

abstract type AbstractCrystallineBoule{T <: SSDFloat} end

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

get_unitful_property(boule::CrystallineBoule, prop::Symbol) = get_unitful_property(boule, Val(prop))

get_unitful_property(boule::CrystallineBoule, ::Val{:impurity_model_parameters}) = boule.impurity_model_parameters .* fit_parameter_units(nameof(boule.impurity_model))

get_impurity_density(boule::CrystallineBoule, z::Number) = boule.impurity_model(to_internal_units(z), boule.impurity_model_parameters) * internal_impurity_quantity

#need to be consistent with get_ functions (do they always return unitfill values?)
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
# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

"""
    AbstractDetectorDesign{T <: SSDFloat, G <: AbstractDesignGeometry}

Supertype for detector designs — a geometry plus the simulated characterization
results (depletion voltage, operating voltage, minimum bulk field, mass).
"""
abstract type AbstractDetectorDesign{T <: SSDFloat, G <: AbstractDesignGeometry} end

"""
    DetectorDesign{T,G} <: AbstractDetectorDesign{T,G}

A geometry-typed detector design plus the slots that [`characterize!`](@ref)
fills in after running the field solver. Construct via the convenience
constructors such as [`InvertedCoaxDesign`](@ref).

# Fields
- `name::AbstractString`: detector identifier.
- `geometry::AbstractDesignGeometry`: physical geometry, with type parameter
  `G` carrying the geometry's validity tag.
- `offset::T`: axial position of the p⁺ contact relative to the boule seed
  end, in mm.
- `mass::Union{T, Missing}`: physical crystal mass in g, derived from the
  geometry volume and [`ge_76_density`](@ref). `missing` for invalid geometries.
- `is_simulated::Bool`: whether a field-solver run has populated the
  characterization results below.
- `Emin::Union{T, Missing}`, `Emin_pos::Union{Tuple{T,T}, Missing}`: minimum
  bulk electric field and its `(r, z)` location, in V/cm and mm.
- `Vdep::Union{T, Missing}`: estimated depletion voltage in V.
- `Vop::Union{T, Missing}`: operational voltage used during characterization, V.
"""
mutable struct DetectorDesign{T,G} <: AbstractDetectorDesign{T,G}
    name::AbstractString
    geometry::AbstractDesignGeometry{T}
    offset::T
    mass::Union{T, Missing}
    is_simulated::Bool
    Emin::Union{T, Missing}
    Emin_pos::Union{Tuple{T,T}, Missing}
    Vdep::Union{T, Missing}
    Vop::Union{T, Missing}
end

"""
    DetectorDesign(det::DetectorDesign{T}, geo::AbstractDesignGeometry{T}) -> DetectorDesign

Rebuild `det` with a new geometry `geo`, preserving the `name` and `offset`,
recomputing `mass` from `geo`, and **resetting** the simulation flags
(`is_simulated`, `Emin`, `Vdep`, `Vop`) — the previous characterization is no
longer valid for a different geometry.
"""
function DetectorDesign(det::DetectorDesign{T}, geo::AbstractDesignGeometry{T}) where {T <: SSDFloat}
    DetectorDesign{T, typeof(geo)}(det.name, geo, det.offset, T(ge_76_density)*get_physical_volume(geo), false, missing, missing, missing, missing)
end

"""
    InvertedCoaxDesign(det::DetectorDesign{T}, geo::InvertedCoaxGeometry{T})

ICPC-specific alias for the geometry-swapping [`DetectorDesign(det, geo)`](@ref).
"""
InvertedCoaxDesign(det::DetectorDesign{T}, geo::InvertedCoaxGeometry{T}) where {T <: SSDFloat} = DetectorDesign(det, geo)

"""
    InvertedCoaxDesign(::Type{T}; name, height, radius, pc_radius,
        groove_depth = 2, groove_inner_radius = pc_radius, groove_outer_radius = groove_inner_radius + 3,
        borehole_pc_gap, borehole_radius,
        borehole_taper_height = 0, borehole_taper_angle = 5,
        top_taper_height = 1, top_taper_angle = 45,
        bottom_taper_height = 1, bottom_taper_angle = 45,
        offset, dead_layer_depth = 1) -> DetectorDesign

High-level constructor: build an [`InvertedCoaxGeometry`](@ref) from the named
dimensions and wrap it in a [`DetectorDesign`](@ref) of element type `T`.

All dimensions accept unitful or bare numbers; bare numbers are interpreted in
the package's internal units (mm for lengths, degrees for angles). The
returned design is tagged `ValidGeometry` only when the geometry passes its
internal validity checks **and** `offset >= height` (i.e., the detector fits
inside the boule); otherwise it is tagged `InvalidGeometry` and `mass` is
`missing`.

Defaults reflect typical LEGEND ICPC values: 2 mm groove depth, 3 mm groove
width past the point contact, 5° borehole taper, 45° outer taper at the
corners, 1 mm dead layer. `offset` is the position of the p⁺ point contact
from the seed end of the boule.
"""
function InvertedCoaxDesign(::Type{T};
        name::AbstractString,
        height::Number,
        radius::Number,
        pc_radius::Number,
        groove_depth::Number = 2,
        groove_inner_radius::Number = to_internal_length_units(pc_radius),
        groove_outer_radius::Number = to_internal_length_units(groove_inner_radius) + 3,
        borehole_pc_gap::Number,
        borehole_radius::Number,
        borehole_taper_height::Number = 0,
        borehole_taper_angle::Number = 5,
        top_taper_height::Number = 1,
        top_taper_angle::Number = 45,
        bottom_taper_height::Number = 1,
        bottom_taper_angle::Number = 45,
        offset::Number,
        dead_layer_depth::Number = 1,
    ) where {T <: SSDFloat}

    height = T(to_internal_length_units(height))
    radius = T(to_internal_length_units(radius))
    pc_radius = T(to_internal_length_units(pc_radius))
    groove_depth = T(to_internal_length_units(groove_depth))
    groove_inner_radius = T(to_internal_length_units(groove_inner_radius))
    groove_outer_radius = T(to_internal_length_units(groove_outer_radius))
    borehole_pc_gap = T(to_internal_length_units(borehole_pc_gap))
    borehole_radius = T(to_internal_length_units(borehole_radius))
    borehole_taper_height = T(to_internal_length_units(borehole_taper_height))
    top_taper_height = T(to_internal_length_units(top_taper_height))
    bottom_taper_height = T(to_internal_length_units(bottom_taper_height))
    offset = T(to_internal_length_units(offset))
    dead_layer_depth = T(to_internal_length_units(dead_layer_depth))

    borehole_taper_angle = T(to_internal_angle_units(borehole_taper_angle))
    top_taper_angle = T(to_internal_angle_units(top_taper_angle))
    bottom_taper_angle = T(to_internal_angle_units(bottom_taper_angle))
 
    geo = InvertedCoaxGeometry{T}(
        height, radius, pc_radius, groove_depth, groove_outer_radius, groove_inner_radius, borehole_pc_gap, borehole_radius, borehole_taper_height, borehole_taper_angle, top_taper_height, top_taper_angle, bottom_taper_height, bottom_taper_angle, dead_layer_depth
    )
    if is_valid_geometry(geo) && offset >= geo.height
        DetectorDesign{T, InvertedCoaxGeometry{T,ValidGeometry}}(name, geo, offset, T(ge_76_density)*get_physical_volume(geo), false, missing, missing, missing, missing)
    else
        DetectorDesign{T, InvertedCoaxGeometry{T,InvalidGeometry}}(name, geo, offset, missing, false, missing, missing, missing, missing)
    end
end

function print(io::IO, det::DetectorDesign{T}) where {T <: SSDFloat}
    g1,g2,g3,g4,g5,g6,g7 = get_unicode_rep(det.geometry)
    vopstr = ismissing(det.Vop) ? missing : @sprintf("%.1f V", det.Vop)
    vstr = ismissing(det.Vdep) ? missing : @sprintf("%.1f V", det.Vdep)
    o,vstr = det.is_simulated && ismissing(det.Vdep) ? ("◯", "> $(det.Vop) V") : (" " , vstr)
    estr = ismissing(det.Emin) ? missing : @sprintf("%.1f V/cm @ r = %.1f, z = %.1f", det.Emin, det.Emin_pos...)
    println(io, "$g1  DetectorDesign{$T} - $(det.name)")
    println(io, "$g2  ╰─Geometry: $(typeof(det.geometry))")
    println(io, "$g3  ╰─Offset of p⁺contact from seed end: $(round(1.0*det.offset, digits = 1)*internal_length_unit)")
    println(io, "$g4  ╰─Minimum E field: $estr")
    println(io, "$g5  $o  $g6  ╰─Depletion/Operational Voltage: $vstr / $vopstr")
    println(io, "$g7  ╰─Mass: $( ismissing(det.mass) ? missing : Int(round(det.mass, digits = 0))*internal_mass_unit)")
end

function print(io::IO, det::DetectorDesign{T, <:AbstractDesignGeometry{T,InvalidGeometry}}) where {T <: SSDFloat}
   println(io, "INVALID GEOMETRY!") 
end

function show(io::IO, det::DetectorDesign)
    print(io, det)
end

function show(io::IO, ::MIME"text/plain", det::DetectorDesign)
    show(io, det)
end

"""
    design_to_meta(det::DetectorDesign{T}) -> PropDict

Serialize `det` to the LEGEND detector-metadata `PropDict` layout. Delegates
to the geometry-level [`geo_to_meta`](@ref), forwarding `det.Vop` and
`det.name`.
"""
function design_to_meta(det::DetectorDesign{T})::PropDict where {T}
    geo_to_meta(det.geometry, Vop = det.Vop, name = det.name)
end

"""
    meta_to_design(::Type{T}, meta::PropDict) -> DetectorDesign{T}

Inverse of [`design_to_meta`](@ref): build a [`DetectorDesign`](@ref) of
element type `T` from LEGEND detector metadata. Currently supports
`meta.type == "icpc"` only; throws `ArgumentError` otherwise.

The result is **not yet simulated** (`is_simulated = false`, `Emin = missing`)
but carries any `depletion_voltage_in_V` / `recommended_voltage_in_V` found in
`meta.characterization.manufacturer` as `Vdep` / `Vop`. `offset` defaults to
0 — set it explicitly if you need the boule-cutting position.
"""
function meta_to_design(::Type{T}, meta::PropDict) where {T}
    if meta.type == "icpc"
        geo = meta_to_geo(InvertedCoaxGeometry{T}, meta)
        Vdep = :depletion_voltage_in_V in keys(meta.characterization.manufacturer) ? T(meta.characterization.manufacturer.depletion_voltage_in_V) : missing
        Vop = :recommended_voltage_in_V in keys(meta.characterization.manufacturer) ? T(meta.characterization.manufacturer.recommended_voltage_in_V) : missing
        DetectorDesign{T, typeof(geo)}(meta.name, geo, T(0), T(ge_76_density)*get_physical_volume(geo), false, missing, missing, Vdep, Vop)
    else
        throw(ArgumentError("Unknown geometry type $(meta.type)"))
    end
end

"""
    get_unitful_property(det::DetectorDesign, prop::Symbol)

Return a field of `det` re-attached with its physical units, via `Val`
dispatch on `prop`. Currently supports `:offset` (returns mm).
"""
get_unitful_property(det::DetectorDesign, prop::Symbol) = get_unitful_property(det, Val(prop))

get_unitful_property(det::DetectorDesign, ::Val{:offset}) = det.offset * internal_length_unit

"""
    SolidStateDetectors.Simulation{T}(det::DetectorDesign{T}, boule::CrystallineBoule{T},
                                       env::HPGeEnvironment = HPGeEnvironment(); kwargs...)
    SolidStateDetectors.Simulation{T}(det::DetectorDesign{T}, imp_model::AbstractImpurityDensity{T},
                                       env::HPGeEnvironment = HPGeEnvironment(); kwargs...)
    SolidStateDetectors.SolidStateDetector{T}(det::DetectorDesign{T}, imp_model::AbstractImpurityDensity{T},
                                              env::HPGeEnvironment = HPGeEnvironment(); kwargs...)
    SolidStateDetectors.SolidStateDetector{T}(det::DetectorDesign{T},
                                              env::HPGeEnvironment = HPGeEnvironment(); kwargs...)

Convert a [`DetectorDesign`](@ref) into the corresponding SolidStateDetectors
object, ready for field-solver use. Internally:

1. [`design_to_meta`](@ref) serializes the geometry to a `PropDict`.
2. SSD builds the simulation / detector from that metadata plus
   [`get_default_xtal_meta`](@ref) as a placeholder impurity profile.
3. The placeholder is replaced via `SolidStateDetector(ssd, imp_model)`:
   - given an `AbstractImpurityDensity`, it is used directly;
   - given a [`CrystallineBoule`](@ref), the impurity model is built on the
     fly from `boule.impurity_model` and `boule.impurity_model_parameters`,
     evaluated at the design's cutting offset (`det.offset`) — convenience
     for the common case where the boule was already fit.

The `SolidStateDetector` method without an impurity argument leaves the
placeholder constant in place and is useful for geometry-only inspection.

`kwargs` are forwarded to SSD's underlying constructor (e.g. `T` precision
overrides, `verbose`, etc.).
"""
function SolidStateDetectors.Simulation{T}(det::DetectorDesign{T}, imp_model::AbstractImpurityDensity{T}, env::HPGeEnvironment = HPGeEnvironment(); kwargs...) where {T<:AbstractFloat}
    meta = design_to_meta(det)
    sim = Simulation{T}(LegendData, meta, get_default_xtal_meta(det), env, verbose = false; kwargs...)
    sim.detector = SolidStateDetector(sim.detector, imp_model)
    sim
end

function SolidStateDetectors.SolidStateDetector{T}(det::DetectorDesign{T}, imp_model::AbstractImpurityDensity{T}, env::HPGeEnvironment = HPGeEnvironment(); kwargs...) where {T<:AbstractFloat}
    meta = design_to_meta(det)
    ssd = SolidStateDetector{T}(LegendData, meta, get_default_xtal_meta(det), env, verbose = false; kwargs...)
    SolidStateDetector(ssd, imp_model)
end

function SolidStateDetectors.SolidStateDetector{T}(det::DetectorDesign{T}, env::HPGeEnvironment = HPGeEnvironment(); kwargs...) where {T<:AbstractFloat}
    meta = design_to_meta(det)
    SolidStateDetector{T}(LegendData, meta, get_default_xtal_meta(det), env, verbose = false; kwargs...)
end
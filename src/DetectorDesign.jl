# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

abstract type AbstractDetectorDesign{T <: SSDFloat, G <: AbstractDesignGeometry} end

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
    offset::Number
    ) where {T}

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

    borehole_taper_angle = T(to_internal_angle_units(borehole_taper_angle))
    top_taper_angle = T(to_internal_angle_units(top_taper_angle))
    bottom_taper_angle = T(to_internal_angle_units(bottom_taper_angle))

    rpass =  pc_radius ≤ groove_inner_radius < groove_outer_radius < radius && borehole_radius < radius
    hpass = borehole_taper_height ≤ height - borehole_pc_gap && top_taper_height < height 
    V = if rpass && hpass
        ValidGeometry
    else
        InvalidGeometry
    end
    geo = InvertedCoaxGeometry{T, V}(
        height, radius, pc_radius, groove_depth, groove_outer_radius,
        groove_inner_radius, borehole_pc_gap, borehole_radius, borehole_taper_height,
        borehole_taper_angle, top_taper_height, top_taper_angle, bottom_taper_height, bottom_taper_angle
    )
    DetectorDesign{T, typeof(geo)}(name, geo, offset, T(ge_76_density)*get_physical_volume(geo), false, missing, missing, missing, missing)
end

function print(io::IO, det::DetectorDesign{T}) where {T <: SSDFloat}
    g1,g2,g3,g4,g5,g6,g7 = get_unicode_rep(det.geometry)
    vstr = ismissing(det.Vdep) ? missing : Int(round(det.Vdep, digits = 0))*internal_voltage_unit 
    o,vstr = det.is_simulated && ismissing(det.Vdep) ? ("◯", "> $(det.Vop) V") : (" " , vstr)
    estr = ismissing(det.Emin) ? missing : "$(Int(round(det.Emin, digits = 0))*internal_efield_unit) @ r = $(round(1.0*det.Emin_pos[1], digits = 1)*internal_length_unit), z = $(round(1.0*det.Emin_pos[2], digits = 1)*internal_length_unit)"
    println(io, "$g1  DetectorDesign{$T} - $(det.name)")
    println(io, "$g2  ╰─Geometry: $(typeof(det.geometry))")
    println(io, "$g3  ╰─Offset of p⁺contact from seed end: $(round(1.0*det.offset, digits = 1)*internal_length_unit)")
    println(io, "$g4  ╰─Minimum E field: $estr")
    println(io, "$g5  $o  $g6  ╰─Depletion Voltage: $vstr")
    println(io, "$g7  ╰─Mass: $( ismissing(det.mass) ? missing : Int(round(det.mass, digits = 0))*internal_mass_unit)")
end

function print(io::IO, det::DetectorDesign{T, <:AbstractDesignGeometry{T,InvalidGeometry}}) where {T <: SSDFloat}
   println(io, "INVALID GEOMETRY!") 
end

function show(io::IO, det::DetectorDesign{T}) where {T <: SSDFloat}
    print(io, det)
end

function show(io::IO, ::MIME"text/plain", det::DetectorDesign{T}) where {T <: SSDFloat}
    show(io, det)
end

function design_2_meta(det::DetectorDesign{T})::PropDict where {T}
    design_2_meta(det.geometry, Vop = det.Vop, name = det.name)   
end

function SolidStateDetectors.Simulation{T}(det::DetectorDesign{T}, imp_model::AbstractImpurityDensity{T}, env::HPGeEnvironment = HPGeEnvironment()) where {T<:AbstractFloat}
    meta = design_2_meta(det)
    sim = Simulation{T}(LegendData, meta, PropDict(), env)
    sim.detector = SolidStateDetector(sim.detector, imp_model);
    sim
end
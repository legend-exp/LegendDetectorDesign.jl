# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).
struct InvertedCoaxGeometry{T,V} <: AbstractDesignGeometry{T,V}
    height::T
	radius::T
	pc_radius::T
    groove_depth::T
	groove_outer_radius::T
    groove_inner_radius::T
    borehole_pc_gap::T
    borehole_radius::T
	borehole_taper_height::T
	borehole_taper_angle::T
    top_taper_height::T
	top_taper_angle::T
    bottom_taper_height::T
	bottom_taper_angle::T
end

function InvertedCoaxGeometry{T}(
        height::T,
        radius::T,
        pc_radius::T,
        groove_depth::T,
        groove_outer_radius::T,
        groove_inner_radius::T,
        borehole_pc_gap::T,
        borehole_radius::T,
        borehole_taper_height::T,
        borehole_taper_angle::T,
        top_taper_height::T,
        top_taper_angle::T,
        bottom_taper_height::T,
        bottom_taper_angle::T
    ) where {T}

    rpass =  pc_radius ≤ groove_inner_radius < groove_outer_radius < radius && borehole_radius < radius
    hpass = borehole_taper_height ≤ height - borehole_pc_gap && top_taper_height < height

    InvertedCoaxGeometry{T, rpass && hpass ? ValidGeometry : InvalidGeometry}(
        height, radius, pc_radius, groove_depth, groove_outer_radius,
        groove_inner_radius, borehole_pc_gap, borehole_radius, borehole_taper_height,
        borehole_taper_angle, top_taper_height, top_taper_angle, bottom_taper_height, bottom_taper_angle
    )
end

function InvertedCoaxGeometry(geo::InvertedCoaxGeometry{T}; 
        height::T = geo.height,
        radius::T = geo.radius,
        pc_radius::T = geo.pc_radius,
        groove_depth::T = geo.groove_depth,
        groove_outer_radius::T = geo.groove_outer_radius,
        groove_inner_radius::T = geo.groove_inner_radius,
        borehole_pc_gap::T = geo.borehole_pc_gap,
        borehole_radius::T = geo.borehole_radius,
        borehole_taper_height::T = geo.borehole_taper_height,
        borehole_taper_angle::T = geo.borehole_taper_angle,
        top_taper_height::T = geo.top_taper_height,
        top_taper_angle::T = geo.top_taper_angle,
        bottom_taper_height::T = geo.bottom_taper_height,
        bottom_taper_angle::T = geo.bottom_taper_angle
    ) where {T <: SSDFloat}
    
    InvertedCoaxGeometry{T}(
        height, radius, pc_radius, groove_depth, groove_outer_radius, groove_inner_radius, borehole_pc_gap, borehole_radius, borehole_taper_height, borehole_taper_angle, top_taper_height, top_taper_angle, bottom_taper_height, bottom_taper_angle
    )
end

function get_unicode_rep(geo::InvertedCoaxGeometry{T, ValidGeometry}) where {T}
    "╭───╮ ╭───╮", 
    "│   │ │   │",
    "│   │ │   │",
    "│   ╰─╯   │",
    "│  ",  "  │",
    "╰── ─── ──╯"
end

function get_physical_volume(geo::InvertedCoaxGeometry{T, ValidGeometry})::T where {T}
	V_cyl = π*geo.radius^2 * (geo.height - geo.top_taper_height - geo.bottom_taper_height)

	# volume of the upper part, where there might be outer taper 
	# -> detector there is a truncated cone
	rtop = geo.radius - geo.top_taper_height*tand(geo.top_taper_angle)
	V_top_cone = 1/3 * π * (geo.radius^2+geo.radius*rtop+rtop^2) * geo.top_taper_height 
    
    rbottom = geo.radius - geo.bottom_taper_height*tand(geo.bottom_taper_angle)
	V_bottom_cone = 1/3 * π * (geo.radius^2+geo.radius*rbottom+rbottom^2) * geo.bottom_taper_height

	# volume of the cylindrical part of the borehole
	V_borehole_cyl = π*geo.borehole_radius^2 * ((geo.height - geo.borehole_pc_gap) - geo.borehole_taper_height)

	# volume of the inner tapered part of the borehole 
	# -> a truncated cone
	rtop_borehole = geo.borehole_radius + geo.borehole_taper_height*tand(geo.borehole_taper_angle)
	V_borehole_cone = 1/3 * π * (rtop_borehole^2+rtop_borehole*geo.borehole_radius+geo.borehole_radius^2) * geo.borehole_taper_height

	# total volume 
	Vtot = V_cyl + V_top_cone + V_bottom_cone - V_borehole_cyl - V_borehole_cone
    ustrip(u"cm^3", Vtot*internal_length_unit^3)
end

function design_2_meta(geo::InvertedCoaxGeometry{T,ValidGeometry}; 
        Vop::Union{Missing, <:Real} = missing, 
        name::AbstractString = "IC"
    )::PropDict where {T}

    geodict = PropDict(
        "height_in_mm" => geo.height,
        "radius_in_mm" => geo.radius,
        "borehole" => PropDict(
            "radius_in_mm" => geo.borehole_radius, 
            "depth_in_mm" => geo.height - geo.borehole_pc_gap),
        "groove" => PropDict(
            "depth_in_mm" => geo.groove_depth, 
            "radius_in_mm" => PropDict("outer" => geo.groove_outer_radius, "inner" => geo.groove_inner_radius)),
        "pp_contact" => PropDict("radius_in_mm" => geo.pc_radius, "depth_in_mm" => 0),
        "taper" => PropDict(
            "top" => PropDict("angle_in_deg" => geo.top_taper_angle, "height_in_mm" => geo.top_taper_height),
            "bottom" => PropDict("angle_in_deg" => geo.bottom_taper_angle, "height_in_mm" => geo.bottom_taper_height),
            "borehole" => PropDict("angle_in_deg" => geo.borehole_taper_angle, "height_in_mm" => geo.borehole_taper_height))
    ) 
    PropDict(
        "geometry" => geodict,
        "type" => "icpc",
        "name" => name,
        "characterization" => PropDict("manufacturer" => PropDict("recommended_voltage_in_V" => ismissing(Vop) ? default_operational_V : Vop))
    ) 
end

function print(io::IO, geo::InvertedCoaxGeometry{T, ValidGeometry}) where {T <: SSDFloat}
    g1,g2,g3,g4,g5,g6,g7 = get_unicode_rep(geo)
    br, pr,r, h, hg = Int.(round.((geo.borehole_radius, geo.pc_radius, geo.radius, geo.height, geo.borehole_pc_gap)))
    println(io, "     ╭╮$(lpad(string(br), 2, ' '))ₘₘ       $(typeof(geo))")
    println(io, "$g1  ╮    ╰─Bulk: Radius, Height, p⁺contact radius:")
    println(io, "$g2  │      ╰─$(geo.radius), $(geo.height), $(geo.pc_radius) mm")
    println(io, "$g3  │    ╰─Borehole: Gap, Radius, Taper height/angle:")
    println(io, "$g4  │      ╰─$(geo.borehole_pc_gap), $(geo.borehole_radius), $(geo.borehole_taper_height) mm / $(geo.borehole_taper_angle)°")
    println(io, "$g5$(lpad(string(hg), 3, ' '))ₘₘ$g6 $(lpad(string(h), 3, ' '))ₘₘ ╰─Goove: Depth, Inner/Outer radius:")
    println(io, "$g7  ╯      ╰─$(geo.groove_depth), $(geo.groove_inner_radius) / $(geo.groove_outer_radius) mm")
    println(io, "     ╰─╯ $(lpad(string(pr), 2, ' '))ₘₘ     ╰─Outer Taper: Top, Bottom height/angle:")
    println(io, "     ╰────╯ $(lpad(string(r), 2, ' '))ₘₘ    ╰─$(geo.top_taper_height) mm / $(geo.top_taper_angle)°, $(geo.bottom_taper_height) mm / $(geo.bottom_taper_angle)°")
    #println(io, "     ╰────╯ $(lpad(string(r), 2, ' '))ₘₘ")
end


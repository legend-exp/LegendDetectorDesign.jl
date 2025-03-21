# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

module LegendDetectorDesignPlotsExt

using RecipesBase
import Plots
using Printf
using LinearAlgebra
using LegendDetectorDesign
using SolidStateDetectors

abstract type AbstractMeasurement{T <: Real} end

struct LinearMeasurement{T} <: AbstractMeasurement{T}
    p1::Tuple{T,T}
    p2::Tuple{T,T}
    offset::T
    annotation::AbstractString
    guides::Tuple{T,T}
    guide_offset::T
end

function LinearMeasurement{T}(p1::Tuple{<:Real,<:Real}, p2::Tuple{<:Real,<:Real}, offset::Real, annotation::AbstractString) where {T}
    LinearMeasurement{T}(p1, p2, offset, annotation, (offset,offset), 1)
end

struct AngularMeasurement{T} <: AbstractMeasurement{T}
    origin::Tuple{T,T}
    reference::T
    offset::T
    annotation::AbstractString
    guides::Tuple{T,T}
end

get_internal_unit(m::LinearMeasurement) = internal_length_unit
get_internal_unit(m::AngularMeasurement) = internal_angle_unit


function get_annotation_info(m::LinearMeasurement{T}) where {T}
    p1 = [m.p1...]
    p2 = [m.p2...]
    u = p2 - p1
    unorm = norm(u)
    u = u/unorm
    p1 = p1 + m.offset*[−u[2],u[1]]
    p2 = p2 + m.offset*[−u[2],u[1]]
    (p1+p2)/2 + 3*[−u[2],u[1]], atand(u[2], u[1]), unorm
end

@recipe function f(m::LinearMeasurement{T}; arrows_and_guides = :both, annotate_length = true, measurementfontsize = 6) where {T}
    seriestype := :shape
    linecolor --> :black
    fillcolor := :black
    aspect_ratio := 1.0
    label --> nothing
    p1 = [m.p1...]
    p2 = [m.p2...]
    u = p2 - p1
    unorm = norm(u)
    u = u/unorm
    θ = atan(u[2], u[1])

    p1 = p1 + m.offset*[−u[2],u[1]]
    p2 = p2 + m.offset*[−u[2],u[1]]

    mid = (p1+p2)/2 + 3*[−u[2],u[1]]

    atext = if annotate_length 
        alabel = m.annotation * (unorm % 1 > 0 ? "$(@sprintf("%.1f", unorm))" : "$(@sprintf("%.0f", unorm))")
        Plots.text(alabel, measurementfontsize, :black, :center, rotation=rad2deg(θ))
    else
        Plots.text(m.annotation, measurementfontsize, :black, :center, rotation=rad2deg(θ))
    end

    annotations --> (mid[1], mid[2], atext)
        
    R = T.([cos(θ) -sin(θ); sin(θ)  cos(θ)])
    a_len = 1
    invert_offset, extend = if unorm < a_len*8 (unorm, 5) else (0,0) end
    l1 = p1 - extend*a_len*u
    l2 = p2 + extend*a_len*u
    line = [l1[1], l2[1]], [l1[2], l2[2]]

    plotelements = [line]
    if arrows_and_guides == :end || arrows_and_guides == :both
        arrow_points = [T.([-a_len,-0.5]), T.([a_len, 0]), T.([-a_len,0.5])]
        arrow_points_t = broadcast(p -> R*p + p2 - a_len*u - invert_offset*u, arrow_points)
        arrow = getindex.(arrow_points_t,1),getindex.(arrow_points_t, 2)
        
        guide_points = [T.([0,sign(m.offset)]), T.([0,-m.guides[2]+m.guide_offset*sign(m.offset)])]
        guide_points_t = broadcast(p -> R*p + p2, guide_points)
        guide = getindex.(guide_points_t,1), getindex.(guide_points_t,2)
        push!(plotelements, arrow, guide)
    end
    if arrows_and_guides == :start || arrows_and_guides == :both
        arrow_points = [T.([a_len,-0.5]), T.([-a_len, 0]), T.([a_len,0.5])]
        arrow_points_t = broadcast(p -> R*p + p1 + a_len*u + invert_offset*u, arrow_points)
        arrow = getindex.(arrow_points_t,1),getindex.(arrow_points_t,2)

        guide_points = [T.([0,sign(m.offset)]), T.([0,-m.guides[1]+m.guide_offset*sign(m.offset)])]
        guide_points_t = broadcast(p -> R*p + p1, guide_points)
        guide = getindex.(guide_points_t,1), getindex.(guide_points_t,2)
        push!(plotelements, arrow, guide)
    end
    plotelements
end

@recipe function f(vm::AbstractVector{<:LinearMeasurement}; annotate_length = true, measurementfontsize = 6)
    aspect_ratio := 1.0
    arrows_and_guides --> :both

    info = get_annotation_info.(vm)
    loc = getindex.(info, 1)
    rot = getindex.(info, 2)
    len = getindex.(info, 3)
    ant = getproperty.(vm, :annotation)

    alen = annotate_length ? broadcast(l -> l % 1 > 0 ? "$(@sprintf("%.1f", l))" : "$(@sprintf("%.0f", l))", len) : ""
    labels = ant .* alen 

    labelsf = [Plots.text(labels[i], measurementfontsize, :black, :center, rotation=rot[i]) for i in 1:length(vm)]
    for m in vm
        @series begin
            m
        end
    end

    @series begin
        seriestype := :scatter
        markersize := 0
        markercolor := :white
        series_annotations := labelsf
        getindex.(loc, 1), getindex.(loc, 2)
    end
end

@recipe function f(geo::InvertedCoaxGeometry{T, ValidGeometry}; corner_rounding = true, include_measurements = false) where {T}
    fc = include_measurements ? :white : :grey
    seriestype := :shape
    linecolor --> :black
    fillcolor --> fc
    aspect_ratio := 1.0
    label --> nothing
    H = geo.height
    R = geo.radius

    if include_measurements
        ticks := false
        guide := ""
        axis := false
    end
    
    bottom_taper_x = [R-geo.bottom_taper_height*tand(geo.bottom_taper_angle),R]
    bottom_taper_y = [0,geo.bottom_taper_height]
    
    if corner_rounding && geo.bottom_taper_angle == 45
        x = range(0, geo.bottom_taper_height, 25)
        bottom_taper_x = x .+ bottom_taper_x[1]
        bottom_taper_y = bottom_taper_y[2] .- sqrt.(bottom_taper_y[2]^2 .- x.^2)
    end
    
    detslicebot_x = [0,geo.groove_inner_radius,geo.groove_inner_radius,geo.groove_outer_radius,geo.groove_outer_radius,bottom_taper_x...]
    detslicebot_y = [0,0,geo.groove_depth,geo.groove_depth,0,bottom_taper_y...]
    detbot_y = vcat(reverse(detslicebot_y), detslicebot_y)
    detbot_x = vcat(-reverse(detslicebot_x), detslicebot_x)
    
    top_taper_x = [R,R-geo.top_taper_height*tand(geo.top_taper_angle)]
    top_taper_y = [H-geo.top_taper_height, H]
    
    if corner_rounding && geo.top_taper_angle == 45
        x = range(geo.top_taper_height, 0, 25)
        top_taper_x = x .+ top_taper_x[2]
        top_taper_y = top_taper_y[1] .+ sqrt.(geo.top_taper_height^2 .- x.^2)
    end
    
    borehole_x = geo.borehole_radius+geo.borehole_taper_height*tand(geo.borehole_taper_angle)
    detslicetop_x = [top_taper_x...,borehole_x,geo.borehole_radius,geo.borehole_radius,geo.borehole_radius]
    detslicetop_y = [top_taper_y...,H,H-geo.borehole_taper_height,geo.borehole_pc_gap, geo.borehole_pc_gap]
    
    dettop_y = vcat(detslicetop_y, reverse(detslicetop_y))
    dettop_x = vcat(detslicetop_x, -reverse(detslicetop_x))
    
    det_y = vcat(detbot_y, dettop_y, detbot_y[1:1])
    det_x = vcat(detbot_x, dettop_x, detbot_x[1:1])
    @series begin
        linestyle := :solid
        det_x, det_y
    end
    if include_measurements
        @series begin
            annotate_length --> true
            linestyle := :solid
            [
                LinearMeasurement{T}((-R, 0), (-R, H),16,""),
                LinearMeasurement{T}((-R, H), (R, H),16,"Ø"),
                LinearMeasurement{T}((-geo.groove_inner_radius, 0), (geo.groove_inner_radius, 0),-10,"Ø"),
                LinearMeasurement{T}((-geo.groove_outer_radius, 0), (geo.groove_outer_radius, 0),-20,"Ø"),
                LinearMeasurement{T}((-geo.borehole_radius, geo.borehole_pc_gap), (geo.borehole_radius,  geo.borehole_pc_gap),-9,"Ø"),
                LinearMeasurement{T}((-borehole_x, geo.borehole_pc_gap), (-borehole_x,  H),10,"",(10+borehole_x-geo.borehole_radius,10), 2),
                LinearMeasurement{T}((borehole_x, H - geo.borehole_taper_height), (borehole_x, H),-10,"",(-10-borehole_x+geo.borehole_radius,-10), 2),
                LinearMeasurement{T}((-borehole_x, H), (borehole_x,  H),8,"Ø"),
                LinearMeasurement{T}((-geo.groove_outer_radius, 0), (-geo.groove_outer_radius, geo.groove_depth),8,"")
            ]
        end
        if corner_rounding
            @series begin
                linestyle := :solid
                annotate_length := false
                arrows_and_guides := :start
                o = 10
                [
                    LinearMeasurement{T}((R, H),(R+o, H + o), 0, "R$(@sprintf("%.0f", geo.top_taper_height))"),
                    LinearMeasurement{T}((R, 0), (R+o, -o), 0, "R$(@sprintf("%.0f", geo.bottom_taper_height))"),
                ]
            end
        else
            ##Outer taper measurements
        end
        @series begin
            linestyle --> :dash
                [
                    [-geo.groove_outer_radius, 0, 1, 2, 3, geo.groove_outer_radius], [-geo.borehole_radius, geo.borehole_radius],
                ],
                [
                    [geo.groove_depth, geo.groove_depth, geo.groove_depth, geo.groove_depth, geo.groove_depth, geo.groove_depth], [H - geo.borehole_taper_height, H - geo.borehole_taper_height]
                ]
            
        end
    end
end

@recipe function f(det::DetectorDesign{T}; crystalname = det.name[4:end-1], technical_drawing = false) where {T}
    size --> (1000,500)
    aspect_ratio := 1.0
    fmt --> :png
    dpi --> 300
    @series begin
        subplot := 1
        corner_rounding --> true
        if technical_drawing 
            include_measurements --> true
            measurementfontsize --> 7
            xlims := (-220, 50)
            ylims := (-30, 135)
        end
        det.geometry
    end
    if technical_drawing 
        @series begin
            xlims := (-220, 50)
            ylims := (-30, 135)
            subplot := 1
            ticks := false
            guide := ""
            axis := false
            seriestype := :scatter
            markersize := 0
            markercolor := :white
            label --> nothing
            series_annotations := [
                Plots.text("Crystal " * crystalname, 10, :black, :left),
                Plots.text("Detector name: " * det.name, 6, :black, :left),
                Plots.text("$(Int(floor(det.mass))) g", 6, :black, :left)
            ]
            [-180, -180, -180], [-18, -23.5, -28]
        end

        detector = SolidStateDetector{T}(det)
        
        @series begin
            subplot := 2
            inset := (1, Plots.bbox(0.0, 0.1, 0.5, 0.9, :bottom, :left))
            ticks := false
            guide := ""
            axis := false
            label --> nothing
            xlims := (-0.05,0.05)
            ylims := (-0.05,0.05)
            zlims := (-0.005,0.11)
            camera --> (20,30)
            seriestype := :samplesurface
            markeralpha --> 0.1
            n_samples --> 80 
            markersize --> 0.5
            markeralpha --> 0.5
            detector.semiconductor
        end
    end
end

end # module

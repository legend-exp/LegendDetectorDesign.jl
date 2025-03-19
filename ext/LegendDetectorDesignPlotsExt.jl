# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

"""
    LegendDetectorDesign

Detector design tools for the LEGEND experiment.
"""
module LegendDetectorDesignPlotsExt

using RecipesBase
import Plots
using Printf
using LinearAlgebra
using LegendDetectorDesign

abstract type AbstractMeausurement{T <: Real} end

struct LinearMeausurement{T} <: AbstractMeausurement{T}
    p1::Tuple{T,T}
    p2::Tuple{T,T}
    offset::T
    annotation::AbstractString
end

function get_annotation_info(m::LinearMeausurement{T}) where {T}
    p1 = [m.p1...]
    p2 = [m.p2...]
    u = p2 - p1
    unorm = norm(u)
    u = u/unorm
    p1 = p1 + m.offset*[−u[2],u[1]]
    p2 = p2 + m.offset*[−u[2],u[1]]
    (p1+p2)/2 + 3*[−u[2],u[1]], atand(u[2], u[1]), unorm
end

@recipe function f(m::LinearMeausurement{T}; at_both_ends = true, guides = sign(m.offset) .* T.((1,-abs(m.offset)+1)), annotate_measurement = true) where {T}
    seriestype := :shape
    linecolor --> :black
    fillcolor := :black
    aspect_ratio := 1.0
    label --> nothing
    
    p1 = [m.p1...]
    p2 = [m.p2...]
    u = p2 - p1
    u = u/norm(u)
    θ = atan(u[2], u[1])

    p1 = p1 + m.offset*[−u[2],u[1]]
    p2 = p2 + m.offset*[−u[2],u[1]]

    mid = (p1+p2)/2 + 3*[−u[2],u[1]]
    atext = Plots.text(m.annotation*"$(@sprintf("%.1f", norm(p2-p1)))", 6, :black, :center, rotation=rad2deg(θ))
    annotations --> (mid[1], mid[2], atext)
        
    R = T.([cos(θ) -sin(θ); sin(θ)  cos(θ)])
    a_len = 1
    arrow_points = [T.([-a_len,-0.5]), T.([a_len, 0]), T.([-a_len,0.5])]
    arrow_points_t = broadcast(p -> R*p + p2 - a_len*u, arrow_points)
    arrow1 = getindex.(arrow_points_t,1),getindex.(arrow_points_t, 2)

    guide_points = [T.([0,guides[1]]), T.([0,guides[2]])]
    guide_points_t = broadcast(p -> R*p + p2, guide_points)
    guide1 = getindex.(guide_points_t,1), getindex.(guide_points_t,2)
    line = [p1[1], p2[1]], [p1[2], p2[2]]
    if at_both_ends
        arrow_points = [T.([a_len,-0.5]), T.([-a_len, 0]), T.([a_len,0.5])]
        arrow_points_t = broadcast(p -> R*p + p1 + a_len*u, arrow_points)
        arrow2 = getindex.(arrow_points_t,1),getindex.(arrow_points_t,2)

        guide_points_t = broadcast(p -> R*p + p1, guide_points)
        guide2 = getindex.(guide_points_t,1), getindex.(guide_points_t,2)
        
        [line, guide1, guide2, arrow1, arrow2]
    else
        [line, guide1, arrow1]
    end
end

@recipe function f(vm::AbstractVector{<:LinearMeausurement})
    aspect_ratio := 1.0

    info = get_annotation_info.(vm)
    loc = getindex.(info, 1)
    rot = getindex.(info, 2)
    len = getindex.(info, 3)
    ant = getproperty.(vm, :annotation)
    labels = [Plots.text(ant[i]*"$(@sprintf("%.1f", len[i]))", 5, :black, :center, rotation=rot[i]) for i in 1:length(vm)]
    for m in vm
        @series begin
            #label := nothing
            m
        end
    end
    @series begin
        seriestype := :scatter
        markersize := 0
        markercolor := :white
        series_annotations := labels
        getindex.(loc, 1), getindex.(loc, 2)
    end
end

@recipe function f(geo::InvertedCoaxGeometry{T, ValidGeometry}; corner_rounding = true, include_measurements = false) where {T}
    fc = include_measurements ? :white : :grey
    seriestype --> :shape
    linecolor --> :black
    fillcolor --> fc
    aspect_ratio --> 1.0
    label --> nothing
    H = geo.height
    R = geo.radius

    if include_measurements
        ticks --> false
        guide --> ""
        axis --> false
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

    #plot(det_x, det_y, color = :black, lw = 1, aspect_ratio = 1, legend = false, ticks = false,
      #guide = "", axis = false)
    @series begin
        det_x, det_y
    end
    if include_measurements
        @series begin
            [
                LinearMeausurement{T}((-geo.radius, 0), (-geo.radius, geo.height),16,""),
                LinearMeausurement{T}((-geo.radius, geo.height), (geo.radius, geo.height),16,"ø"),
                LinearMeausurement{T}((-geo.groove_inner_radius, 0), (geo.groove_inner_radius, 0),-10,"ø"),
                LinearMeausurement{T}((-geo.groove_outer_radius, 0), (geo.groove_outer_radius, 0),-20,"ø"),
                LinearMeausurement{T}((-geo.borehole_radius, geo.borehole_pc_gap), (geo.borehole_radius,  geo.borehole_pc_gap),-9,"ø"),
                LinearMeausurement{T}((-borehole_x, geo.borehole_pc_gap), (-borehole_x,  geo.height),10,""),
                LinearMeausurement{T}((borehole_x, geo.height - geo.borehole_taper_height), (borehole_x,  geo.height),-10,""),
                LinearMeausurement{T}((-borehole_x, geo.height), (borehole_x,  geo.height),8,"ø")
            ]
        end
        #=
        @series begin
            # dotted lines
        end
        =#
    end
end

end # module

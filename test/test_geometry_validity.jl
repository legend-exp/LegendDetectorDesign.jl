# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

using Test
using LegendDetectorDesign

@testset "Geometry validity" begin
    T = Float32
    ldet = InvertedCoaxDesign(T,
        name = "TestDesign",
        height = 105.0,
        radius = 45,
        pc_radius = 12,
        borehole_pc_gap = 30,
        borehole_taper_height = 46,
        borehole_radius = 4,
        offset = 0
        )
    @test typeof(ldet.geometry).parameters[2] == ValidGeometry

    ldet2 = InvertedCoaxDesign(T,
        name = "TestDesign",
        height = 105.0,
        radius = 45,
        pc_radius = 12,
        borehole_pc_gap = 30,
        borehole_taper_height = 86,
        borehole_radius = 4,
        offset = 0
        )
    @test typeof(ldet2.geometry).parameters[2] == InvalidGeometry
end # testset

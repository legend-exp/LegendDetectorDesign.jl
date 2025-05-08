# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

import Test

Test.@testset "Package LegendDetectorDesign" begin
    include("test_aqua.jl")

    # include("test_some_source_file.jl")

    include("test_docs.jl")

    include("test_geometry_validity.jl")
end # testset

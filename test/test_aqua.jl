# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

import Test
import Aqua
import LegendDetectorDesign

Test.@testset "Package ambiguities" begin
    Test.@test isempty(Test.detect_ambiguities(LegendDetectorDesign))
end # testset

Test.@testset "Aqua tests" begin
    Aqua.test_all(
        LegendDetectorDesign,
        ambiguities = true
    )
end # testset

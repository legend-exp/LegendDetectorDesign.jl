# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

using Test
using LegendDetectorDesign
import Documenter

Documenter.DocMeta.setdocmeta!(
    LegendDetectorDesign,
    :DocTestSetup,
    :(using LegendDetectorDesign);
    recursive=true,
)
Documenter.doctest(LegendDetectorDesign)

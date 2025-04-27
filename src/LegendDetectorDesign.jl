# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

"""
    LegendDetectorDesign

Detector design tools for the LEGEND experiment.
"""
module LegendDetectorDesign

using LinearAlgebra
using PropDicts
using Unitful
using LegendDataManagement
using SolidStateDetectors
using Printf

import SolidStateDetectors: 
            SSDFloat, AbstractImpurityDensity, Simulation, AbstractCoordinatePoint, update_till_convergence!, mark_bulk_bits!, mark_undep_bits!

import Base: show, print, println, *

include("Units.jl")
include("Geometry/Geometry.jl")
include("DetectorDesign.jl")
include("ImpurityDensities.jl")
include("ElectricField.jl")
include("Characterize.jl")

export DetectorDesign, InvertedCoaxDesign, InvertedCoaxGeometry, LinBouleImpurityDensity, ParBouleImpurityDensity, LinExpBouleImpurityDensity, ParExpBouleImpurityDensity, ValidGeometry, InvalidGeometry, characterize!

end # module

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
using Interpolations
using StaticArrays
using YAML
using OrderedCollections

import SolidStateDetectors: 
        SSDFloat, AbstractImpurityDensity, Simulation, AbstractCoordinatePoint, update_till_convergence!, mark_bulk_bits!, mark_undep_bits!, _adapt_weighting_potential_to_electric_potential_grid!

import Base: show, print, println

include("Units.jl")
include("Geometry/Geometry.jl")
include("DetectorDesign.jl")
include("ImpurityDensities.jl")
include("ElectricField.jl")
include("CrystallineBoule.jl")
include("Characterize.jl")

export DetectorDesign, InvertedCoaxDesign, InvertedCoaxGeometry, BouleGeometry, CrystallineBoule, LinBouleImpurityDensity, ParBouleImpurityDensity, LinExpBouleImpurityDensity, ParExpBouleImpurityDensity, ValidGeometry, InvalidGeometry, characterize!, fit_function, from_internal_units, to_internal_units, get_unitful_property, get_impurity_density, boule_to_meta

end # module

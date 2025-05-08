# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

function is_local_bulk_Efield_minima(sim::Simulation, x₀::CartesianIndex{3}) #move to SSD?
    @assert !ismissing(sim.point_types) "Please calculate the electric potential first using `calculate_electric_potential!(sim)`"
    @assert !ismissing(sim.electric_field) "Please calculate the electric potential first using `calculate_electric_field!(sim)`"
    ϕ = norm.(sim.electric_field.data)
    ϕ₀ = ϕ[x₀]
    indices = CartesianIndices(sim.electric_potential)
    maxI = last(indices)
    minI = first(indices)
    Δ = oneunit(minI)
    neighbours = max(x₀ - Δ, minI):min(x₀ + Δ, maxI)
    is_loc_min = true
    for x in neighbours
        Δx = sum(abs.((x₀ - x).I))
        if Δx == 1
            is_loc_min &= (ϕ[x] >= ϕ₀ && sim.point_types.data[x] & SolidStateDetectors.bulk_bit .> 0)
        end
    end
    is_loc_min
end

function find_minimum_Efield_in_bulk(sim::Simulation; require_local_min = false, verbose = false) #move to SSD?
    @assert !ismissing(sim.point_types) "Please calculate the electric potential first using `calculate_electric_potential!(sim)`"
    @assert !ismissing(sim.electric_field) "Please calculate the electric potential first using `calculate_electric_field!(sim)`"
    bulk_points = findall(sim.point_types.data .& SolidStateDetectors.bulk_bit .> 0)
    ϕ = norm.(sim.electric_field.data)
    Vcm = uconvert(internal_efield_unit, 1*SolidStateDetectors.internal_efield_unit)
    mm = uconvert(internal_length_unit, 1*SolidStateDetectors.internal_length_unit)
    if require_local_min
        idxs = sortperm(ϕ[bulk_points])
        i = findfirst(p -> is_local_bulk_Efield_minima(sim, p), bulk_points[idxs])
        min, pos = if isnothing(i) 
            if verbose @warn "Electric field has no local minima. Returning global min" end
            ϕ[bulk_points[idxs][1]] * Vcm, sim.electric_field.grid[bulk_points[idxs][1]] .* mm
        else
            ϕ[bulk_points[idxs][i]] * Vcm, sim.electric_field.grid[bulk_points[idxs][i]] .* mm
        end
    else
        i = argmin(ϕ[bulk_points])
        ϕ[bulk_points[i]] * Vcm, sim.electric_field.grid[bulk_points[i]] .* mm
    end
end


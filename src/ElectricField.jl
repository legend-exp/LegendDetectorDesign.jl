# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

"""
    is_local_bulk_Efield_minima(sim::Simulation, x₀::CartesianIndex{3}) -> Bool

Check whether the grid point `x₀` is a local minimum of `‖E‖` within the
detector bulk: all 6 nearest neighbours (face-adjacent in `r/φ/z`) must have
field magnitude ≥ the value at `x₀` and be marked as bulk
(`SolidStateDetectors.bulk_bit`).

Requires `sim.point_types` and `sim.electric_field` to have been computed.
"""
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

"""
    find_minimum_Efield_in_bulk(sim::Simulation; require_local_min = false, verbose = false)
        -> (Emin::Quantity, pos::Tuple)

Locate the minimum bulk electric field strength in `sim`.

If `require_local_min = false`, returns the *global* minimum of `‖E‖` over all
bulk-flagged points. If `true`, returns the smallest bulk point that also
satisfies [`is_local_bulk_Efield_minima`](@ref) (i.e., a true local minimum,
not just a low value on a monotonic ramp toward a contact); falls back to the
global minimum if no local minimum exists.

Returns the magnitude in `V/cm` ([`internal_efield_unit`](@ref)) and its
spatial position in the simulation's internal coordinate units.
"""
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


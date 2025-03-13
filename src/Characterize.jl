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
    if require_local_min
        idxs = sortperm(ϕ[bulk_points])
        i = findfirst(p -> is_local_bulk_Efield_minima(sim, p), bulk_points[idxs])
        min, pos = if isnothing(i) 
            if verbose @warn "Electric field has no local minima. Returning global min" end
            ϕ[bulk_points[idxs][1]]*1e-2u"V/cm", sim.electric_field.grid[bulk_points[idxs][1]].* 1e3u"mm"
        else
            ϕ[bulk_points[idxs][i]] * 1e-2u"V/cm", sim.electric_field.grid[bulk_points[idxs][i]].* 1e3u"mm"
        end
    else
        i = argmin(ϕ[bulk_points])
        ϕ[bulk_points[i]] * 1e-2u"V/cm", sim.electric_field.grid[bulk_points[i]].* 1e3u"mm"
    end
end

function characterize!(det::DetectorDesign{T}, imp_model::AbstractImpurityDensity{T}, env::HPGeEnvironment = HPGeEnvironment(); verbose::Bool = false, Vop::Union{Missing, <:Real} = missing) where {T<:AbstractFloat}
    det.Vop = Vop
    sim = Simulation{T}(det, imp_model, env)
    calculate_electric_potential!(sim, 
        refinement_limits = [0.2, 0.1], 
        depletion_handling = true, 
        verbose = verbose)
    det.is_simulated = true
    
    if is_depleted(sim.point_types)
        if verbose println("\nDetector depletes under $(sim.detector.contacts[2].potential) V, refining...\n") end
        calculate_electric_potential!(sim, 
            refinement_limits = [0.05,0.02], 
            depletion_handling = true, 
            verbose = verbose, 
            initialize = false)
        if is_depleted(sim.point_types)
            Vdep = T(to_internal_units(estimate_depletion_voltage(sim, check_for_depletion = false, verbose = verbose)))
            det.Vdep = Vdep
            det.Vop = Vdep + 500
            sim = Simulation{T}(det, imp_model, env)
            if verbose println("\nSimulating at Vdep + 500V...\n") end
            calculate_electric_potential!(sim, 
                refinement_limits = [0.2, 0.1, 0.05, 0.02], 
                depletion_handling = true, 
                verbose = verbose)
            calculate_electric_field!(sim, n_points_in_φ = 2);
            emin, pos = find_minimum_Efield_in_bulk(sim, require_local_min = true)
            det.Emin = T(to_internal_units(emin))
            det.Emin_pos = (T(to_internal_units(pos[1])), T(to_internal_units(pos[3])))
        end
    end
    sim
end
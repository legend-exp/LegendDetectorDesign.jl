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

function initialize_design!(det::DetectorDesign)
    det.Emin = missing
    det.Emin_pos = missing
    det.Vdep = missing
end

function characterize!(det::DetectorDesign{T}, imp_model::AbstractImpurityDensity{T};
        env::HPGeEnvironment = HPGeEnvironment(),
        verbose::Bool = false, 
        Vmax::Real = default_operational_V, 
        refinement_limits::Vector{<:Real} = [0.2, 0.1, 0.05, 0.02]
    ) where {T<:AbstractFloat}

    det.Vop = Vmax
    sim = Simulation{T}(det, imp_model, env)
    calculate_electric_potential!(sim, 
        refinement_limits = refinement_limits[1:min(2,end)], 
        depletion_handling = true, 
        verbose = verbose)
    det.is_simulated = true
    
    if is_depleted(sim.point_types)
        if verbose @info "Detector depletes under $(sim.detector.contacts[2].potential) V" end
        if length(refinement_limits) > 2
            if verbose @info "Refining Grid ...." end
            calculate_electric_potential!(sim, 
                refinement_limits = refinement_limits[3:end], 
                depletion_handling = true, 
                verbose = verbose, 
                initialize = false)
        end
        if length(refinement_limits) <= 2 || is_depleted(sim.point_types)
            Vdep = T(to_internal_units(estimate_depletion_voltage(sim, check_for_depletion = false, verbose = verbose)))
            det.Vdep = Vdep
            if verbose @info "Simulating at Vdep + 500V...\n" end
            sim = characterize!(det, imp_model, Vdep + 500, env = env, check_for_depletion = false, verbose = verbose, refinement_limits = refinement_limits)
        end
    else
        initialize_design!(det)
    end
    sim
end

function characterize!(det::DetectorDesign{T}, imp_model::AbstractImpurityDensity{T}, Vop::Real;
    env::HPGeEnvironment = HPGeEnvironment(),
    check_for_depletion::Bool = true,
    verbose::Bool = false, 
    refinement_limits::Vector{<:Real} = [0.2, 0.1, 0.05, 0.02]
) where {T<:AbstractFloat}
    if !ismissing(det.Vdep) && Vop < det.Vdep
        initialize_design!(det)
        throw(AssertionError("Operational voltage must be greater than depletion voltage"))
    end
    det.Vop = Vop
    sim = Simulation{T}(det, imp_model, env)
    calculate_electric_potential!(sim, 
        refinement_limits = refinement_limits, 
        depletion_handling = true, 
        verbose = verbose)
    det.is_simulated = true
    if !check_for_depletion || is_depleted(sim.point_types)
        if ismissing(det.Vdep)
            dep_tol = 1e-1
            Vdep = T(to_internal_units(estimate_depletion_voltage(sim, check_for_depletion = false, tolerance = dep_tol, verbose = verbose)))
            if Vop - Vdep <= dep_tol @warn "Detector is not depleted. Depletion voltage is likely incorrect!" end
            det.Vdep = Vdep
        end
        calculate_electric_field!(sim, n_points_in_φ = 2);
        if verbose @info "Looking for minimum electric field in the bulk" end
        emin, pos = find_minimum_Efield_in_bulk(sim, require_local_min = true)
        det.Emin = T(to_internal_units(emin))
        det.Emin_pos = (T(to_internal_units(pos[1])), T(to_internal_units(pos[3])))
        if verbose @info "$(Int(round(det.Emin, digits = 0))*internal_efield_unit) @ r = $(round(1.0*det.Emin_pos[1], digits = 1)*internal_length_unit), z = $(round(1.0*det.Emin_pos[2], digits = 1)*internal_length_unit)" end
    else
        initialize_design!(det)
    end
    sim
end

function characterize!(det::DetectorDesign{T}, imp_model::AbstractImpurityDensity{T}, reference_simulation::Simulation{T}; 
    Vop::Real = reference_simulation.detector.contacts[2].potential,
    verbose::Bool = false
) where {T<:AbstractFloat}

    det.Vop = Vop
    initialize_design!(det)
    reference_simulation.detector = SolidStateDetector{T}(det, imp_model)
    
    update_till_convergence!( reference_simulation, ElectricPotential,
                                    n_iterations_between_checks = 1000,
                                    max_n_iterations = 50000,
                                    depletion_handling = true,
                                    verbose = verbose
                            )
    mark_bulk_bits!(reference_simulation.point_types.data)
    mark_undep_bits!(reference_simulation.point_types.data, reference_simulation.imp_scale.data)
    det.is_simulated = true

    if is_depleted(reference_simulation.point_types)
        Vdep = T(to_internal_units(estimate_depletion_voltage(reference_simulation, check_for_depletion = false, verbose = verbose)))
        det.Vdep = Vdep
        calculate_electric_field!(reference_simulation, n_points_in_φ = 2);
        if verbose @info "Looking for minimum electric field in the bulk" end
        emin, pos = find_minimum_Efield_in_bulk(reference_simulation, require_local_min = true)
        det.Emin = T(to_internal_units(emin))
        det.Emin_pos = (T(to_internal_units(pos[1])), T(to_internal_units(pos[3])))
        if verbose @info "$(Int(round(det.Emin, digits = 0))*internal_efield_unit) @ r = $(round(1.0*det.Emin_pos[1], digits = 1)*internal_length_unit), z = $(round(1.0*det.Emin_pos[2], digits = 1)*internal_length_unit)" end
    else
        initialize_design!(det)
    end
    reference_simulation
end

#### OLD METHODS 
#=
function characterize!(det::DetectorDesign{T}, imp_model::AbstractImpurityDensity{T};
        env::HPGeEnvironment = HPGeEnvironment(),
        verbose::Bool = false, 
        Vop::Real = default_operational_V, 
        refinement_limits::Vector{<:Real} = [0.2, 0.1, 0.05, 0.02]
    ) where {T<:AbstractFloat}
    
    det.Vop = Vop
    det.Emin = missing
    det.Emin_pos = missing
    det.Vdep = missing
    sim = Simulation{T}(det, imp_model, env)
    calculate_electric_potential!(sim, 
        refinement_limits = refinement_limits[1:min(2,end)], 
        depletion_handling = true, 
        verbose = verbose)
    det.is_simulated = true
    
    if is_depleted(sim.point_types)
        if verbose @info "Detector depletes under $(sim.detector.contacts[2].potential) V" end
        if length(refinement_limits) > 2
            if verbose @info "Refining Grid ...." end
            calculate_electric_potential!(sim, 
                refinement_limits = refinement_limits[3:end], 
                depletion_handling = true, 
                verbose = verbose, 
                initialize = false)
        end
        if length(refinement_limits) <= 2 || is_depleted(sim.point_types)
            Vdep = T(to_internal_units(estimate_depletion_voltage(sim, check_for_depletion = false, verbose = verbose)))
            det.Vdep = Vdep
            det.Vop = Vdep + 500
            sim = Simulation{T}(det, imp_model, env)
            if verbose @info "Simulating at Vdep + 500V...\n" end
            calculate_electric_potential!(sim, 
                refinement_limits = refinement_limits, 
                depletion_handling = true, 
                verbose = verbose)
            calculate_electric_field!(sim, n_points_in_φ = 2);
            if verbose @info "Looking for minimum electric field in the bulk" end
            emin, pos = find_minimum_Efield_in_bulk(sim, require_local_min = true)
            det.Emin = T(to_internal_units(emin))
            det.Emin_pos = (T(to_internal_units(pos[1])), T(to_internal_units(pos[3])))
            if verbose @info "$(Int(round(det.Emin, digits = 0))*internal_efield_unit) @ r = $(round(1.0*det.Emin_pos[1], digits = 1)*internal_length_unit), z = $(round(1.0*det.Emin_pos[2], digits = 1)*internal_length_unit)" end
        end
    end
    sim
end
=#
#=
function characterize!(det::DetectorDesign{T}, imp_model::AbstractImpurityDensity{T}, reference_simulation::Simulation{T}; 
    Vop::Real = reference_simulation.detector.contacts[2].potential,
    verbose::Bool = false
) where {T<:AbstractFloat}

    det.Vop = Vop
    det.Emin = missing
    det.Emin_pos = missing
    det.Vdep = missing
    reference_simulation.detector = SolidStateDetector{T}(det, imp_model)
    
    update_till_convergence!( reference_simulation, ElectricPotential,
                                    n_iterations_between_checks = 1000,
                                    max_n_iterations = 50000,
                                    depletion_handling = true,
                                    verbose = verbose
                            )
    mark_bulk_bits!(reference_simulation.point_types.data)
    mark_undep_bits!(reference_simulation.point_types.data, reference_simulation.imp_scale.data)
    
    if is_depleted(reference_simulation.point_types)
        Vdep = T(to_internal_units(estimate_depletion_voltage(reference_simulation, check_for_depletion = false, verbose = verbose)))
        det.Vdep = Vdep
        calculate_electric_field!(reference_simulation, n_points_in_φ = 2);
        if verbose @info "Looking for minimum electric field in the bulk" end
        emin, pos = find_minimum_Efield_in_bulk(reference_simulation, require_local_min = true)
        det.Emin = T(to_internal_units(emin))
        det.Emin_pos = (T(to_internal_units(pos[1])), T(to_internal_units(pos[3])))
        if verbose @info "$(Int(round(det.Emin, digits = 0))*internal_efield_unit) @ r = $(round(1.0*det.Emin_pos[1], digits = 1)*internal_length_unit), z = $(round(1.0*det.Emin_pos[2], digits = 1)*internal_length_unit)" end
    end
    reference_simulation
end
=#
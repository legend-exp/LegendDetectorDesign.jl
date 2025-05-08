# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

function initialize_design!(det::DetectorDesign)
    det.Emin = missing
    det.Emin_pos = missing
    det.Vdep = missing
end

function populate_design!(det::DetectorDesign{T}, sim::Simulation{T}, Vop::T; verbose::Bool = verbose) where {T}
    if ismissing(det.Vdep)
        dep_tol = 1e-1
        Vdep = T(to_internal_units(estimate_depletion_voltage(sim, check_for_depletion = false, tolerance = dep_tol, verbose = verbose)))
        if Vop - Vdep <= dep_tol @warn "Detector is not depleted. Depletion voltage is likely incorrect!" end
        det.Vdep = Vdep
    end
    if ismissing(det.Emin)
        calculate_electric_field!(sim, n_points_in_φ = 2);
        if verbose @info "Looking for minimum electric field in the bulk" end
        emin, pos = find_minimum_Efield_in_bulk(sim, require_local_min = true)
        det.Emin = T(to_internal_units(emin))
        det.Emin_pos = (T(to_internal_units(pos[1])), T(to_internal_units(pos[3])))
        if verbose @info @sprintf("%.2f V/cm @ r = %.2f, z = %.2f", det.Emin, det.Emin_pos...) end
    end
end

function characterize!(det::DetectorDesign{T}, imp_model::AbstractImpurityDensity{T};
        env::HPGeEnvironment = HPGeEnvironment(),
        verbose::Bool = false, 
        Vmax::Real = default_operational_V, 
        refinement_limits::Vector{<:Real} = [0.2, 0.1, 0.05, 0.02]
    ) where {T<:AbstractFloat}

    initialize_design!(det)
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
            Vop = Vdep + 500
            det.Vop = Vop
            sim.detector = SolidStateDetector(sim.detector, contact_id = 2, contact_potential =  Vop)
            calculate_electric_potential!(sim, refinement_limits = refinement_limits, depletion_handling = true, verbose = verbose)
            populate_design!(det, sim, T(Vop), verbose = verbose)
        end
    end
    sim
end

function characterize!(det::DetectorDesign{T}, imp_model::AbstractImpurityDensity{T}, Vop::Real;
    env::HPGeEnvironment = HPGeEnvironment(),
    check_for_depletion::Bool = true,
    verbose::Bool = false, 
    refinement_limits::Vector{<:Real} = [0.2, 0.1, 0.05, 0.02],
    initialize::Bool = true
) where {T<:AbstractFloat}

    if initialize initialize_design!(det) end
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
        populate_design!(det, sim, T(Vop), verbose = verbose)
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
    reference_simulation.detector = SolidStateDetector{T}(det, imp_model) ##is this necessary?
    
    update_till_convergence!( reference_simulation, ElectricPotential,
                                    n_iterations_between_checks = 1000,
                                    max_n_iterations = 50000,
                                    depletion_handling = true,
                                    verbose = verbose
                            )
    mark_bulk_bits!(reference_simulation.point_types.data)
    mark_undep_bits!(reference_simulation.point_types.data, reference_simulation.imp_scale.data)
    det.is_simulated = true

    initialize_design!(det)
    if is_depleted(reference_simulation.point_types)
        populate_design!(det, reference_simulation, T(Vop), verbose = verbose)
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
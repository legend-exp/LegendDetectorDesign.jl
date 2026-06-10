# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

"""
    initialize_design!(det::DetectorDesign)

Clear the simulation-derived fields of `det` (`Emin`, `Emin_pos`, `Vdep`) by
setting them to `missing`. Called at the start of every `characterize!`
pathway that re-solves from scratch, so a partially-populated design from a
previous run never leaks into the next.
"""
function initialize_design!(det::DetectorDesign)
    det.Emin = missing
    det.Emin_pos = missing
    det.Vdep = missing
end

"""
    populate_design!(det::DetectorDesign{T}, sim::Simulation{T}, Vop::T; verbose = false)

Fill the simulation-derived fields of `det` from a converged `sim` at operating
voltage `Vop`:

- `det.Vdep`: estimated depletion voltage via SSD's
  `estimate_depletion_voltage` (only if not already set).
- `det.Emin` and `det.Emin_pos`: minimum bulk field magnitude and its `(r, z)`
  position from [`find_minimum_Efield_in_bulk`](@ref) — `calculate_electric_field!`
  is run with `n_points_in_φ = 2` first.

Emits a warning if `Vop − Vdep` is within the depletion tolerance, i.e. the
detector likely isn't actually depleted at `Vop`.
"""
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

"""
    characterize!(det::DetectorDesign, boule::CrystallineBoule; kwargs...) -> Simulation
    characterize!(det::DetectorDesign, imp_model::AbstractImpurityDensity; Vmax = …, refinement_limits = …, …) -> Simulation
    characterize!(det::DetectorDesign, imp_model::AbstractImpurityDensity, Vop::Real; check_for_depletion = true, initialize = true, …) -> Simulation
    characterize!(det::DetectorDesign, imp_model::AbstractImpurityDensity, reference_simulation::Simulation; Vop = …, verbose = false) -> Simulation

Run a SolidStateDetectors field-solver pass for `det` and populate its
simulation-derived fields (`Vdep`, `Vop`, `Emin`, `Emin_pos`, `mass`,
`is_simulated`). Returns the converged `Simulation`.

The four methods differ in how the impurity density is supplied and how the
solver is driven:

- **`(det, boule)`** — builds the impurity model from the boule's fit
  parameters (offset-adjusted to the detector's cutting position) and
  delegates to the `(det, imp_model)` method. The typical entry point when
  you have a [`CrystallineBoule`](@ref).

- **`(det, imp_model)`** — auto-discover depletion: solves at `Vmax` with two
  refinement levels, then if the detector is depleted, refines further, uses
  the superposition principle (potential + `(Vop − Vbias)·ϕ_W`) to advance to
  `Vdep + 500 V`, and finally populates the design. `refinement_limits` is a
  vector of relative-difference targets; the first two are used for the
  initial sweep and the rest for subsequent refinement.

- **`(det, imp_model, Vop)`** — solve at a fixed bias `Vop`. By default
  resets the design first and only populates results if the detector depletes
  (`check_for_depletion = true`); pass `initialize = false` to retain
  previous `Vdep` / `Vop` and `check_for_depletion = false` to always
  populate.

- **`(det, imp_model, reference_simulation)`** — **warm-start**: reuse the
  grid and converged potential of `reference_simulation`, swap in the new
  detector + impurity model, and re-converge via `update_till_convergence!`
  (no refinement, no cold restart). Much faster when `det` is a small
  perturbation of the reference. Default `Vop` is the reference's current
  contact-2 potential.
"""
function characterize!(det::DetectorDesign{T}, boule::CrystallineBoule{T};
        env::HPGeEnvironment = HPGeEnvironment(),
        verbose::Bool = false,
        Vmax::Real = default_operational_V,
        refinement_limits::Vector{<:Real} = [0.2, 0.1, 0.05, 0.02]
    ) where {T<:AbstractFloat}
    idm = ssd_ptype*impurity_density_model(nameof(boule.impurity_model)){T}(get_unitful_property(boule, :impurity_model_parameters), get_unitful_property(det, :offset))
    characterize!(det, idm, env = env, verbose = verbose, Vmax = Vmax, refinement_limits = refinement_limits)
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
            if verbose @info "Refining Grid ..." end
            calculate_electric_potential!(sim, 
                refinement_limits = refinement_limits[3:end], 
                depletion_handling = true, 
                verbose = verbose, 
                initialize = false)
            if verbose @info "Calculating weighting potential for depletion voltage estimation and projection to Vdep + 500V..." end
            _adapt_weighting_potential_to_electric_potential_grid!(sim, 2)
        end
        if length(refinement_limits) <= 2 || is_depleted(sim.point_types)
            Vdep = T(to_internal_units(estimate_depletion_voltage(sim, check_for_depletion = false, verbose = verbose)))
            det.Vdep = Vdep
            if verbose @info "Using superposition principle to calculate at Vdep + 500V...\n" end
            Vop = Vdep + 500
            det.Vop = Vop
            ϕV = sim.weighting_potentials[2].data
            sim.electric_potential.data .+= (Vop - sim.detector.contacts[2].potential) .* ϕV
            sim.detector = SolidStateDetector(sim.detector, contact_id = 2, contact_potential =  Vop)
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

# Method where a reference simulation is provided, so that the electric potential can be updated till convergence instead of recalculating from scratch at each iteration. The idea is that the det is just a small permulation of the reference simulation, so the electric potential should be close to the reference simulation at each iteration, and thus should converge faster.
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
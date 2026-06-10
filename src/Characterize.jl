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
    adapt_potential_to_depletion_plus_offset!(det::DetectorDesign{T}, sim::Simulation{T};
                                               offset::T = T(500), verbose = false)

Estimate the depletion voltage from `sim`, jump the electric potential to
`Vop = Vdep + offset` analytically via the weighting-potential superposition
trick, and re-tag the bias contact (`contact_id = 2`) to the new `Vop`.

Steps:

1. Adapt the bias contact's weighting potential to the electric-potential
   grid via `_adapt_weighting_potential_to_electric_potential_grid!`.
2. Estimate `Vdep` via SSD's `estimate_depletion_voltage`; store on `det.Vdep`.
3. Use the superposition identity
   `φ(Vop) = φ_solved(V_bias) + (Vop − V_bias) · φ_W`
   to shift the converged electric potential to the new bias — no second
   field solve required. Stores `Vop` on `det.Vop`.
4. Re-tag the detector's contact-2 potential to the new `Vop`.

Does **not** populate `det.Emin` / `det.Emin_pos` — the caller is responsible
for calling [`populate_design!`](@ref) afterwards (typically with `det.Vop`
as the bias).

Shared between the auto-discover-depletion and reference-simulation paths
of [`characterize!`](@ref). `offset` defaults to 500 V to land safely above
the depletion knee.
"""
function adapt_potential_to_depletion_plus_offset!(det::DetectorDesign{T}, sim::Simulation{T}; offset::T = T(500), verbose::Bool = false) where {T}
    if verbose @info "Calculating weighting potential for depletion voltage estimation and projection to Vdep + $(offset) V..." end
    _adapt_weighting_potential_to_electric_potential_grid!(sim, 2)
    Vdep = T(to_internal_units(estimate_depletion_voltage(sim, check_for_depletion = false, verbose = verbose)))
    det.Vdep = Vdep
    if verbose @info "Using superposition principle to calculate at Vdep + $(offset) V...\n" end
    Vop = Vdep + offset
    det.Vop = Vop
    ϕV = sim.weighting_potentials[2].data
    sim.electric_potential.data .+= (Vop - sim.detector.contacts[2].potential) .* ϕV
    sim.detector = SolidStateDetector(sim.detector, contact_id = 2, contact_potential = Vop)
end

"""
    characterize!(det::DetectorDesign, boule::CrystallineBoule, [Vop_or_refsim]; kwargs...) -> Simulation
    characterize!(det::DetectorDesign, imp_model::AbstractImpurityDensity; Vmax = …, refinement_limits = …, …) -> Simulation
    characterize!(det::DetectorDesign, imp_model::AbstractImpurityDensity, Vop::Real; check_for_depletion = true, initialize = true, …) -> Simulation
    characterize!(det::DetectorDesign, imp_model::AbstractImpurityDensity, reference_simulation::Simulation; verbose = false) -> Simulation

Run a SolidStateDetectors field-solver pass for `det` and populate its
simulation-derived fields (`Vdep`, `Vop`, `Emin`, `Emin_pos`, `mass`,
`is_simulated`). Returns the converged `Simulation`.

The methods differ in how the impurity density is supplied and how the
solver is driven:

- **`(det, boule, args...)`** — the typical entry point when you have a
  [`CrystallineBoule`](@ref). Builds the impurity model with
  [`impurity_model_from_boule`](@ref) and forwards `args...` / `kwargs...` to
  whichever `(det, imp_model, args...)` method below matches. So
  `characterize!(det, boule)` runs auto-depletion,
  `characterize!(det, boule, 3000)` solves at a fixed bias, and
  `characterize!(det, boule, sim_ref)` warm-starts from a reference.

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

- **`(det, imp_model, reference_simulation)`** — **warm-start**: deep-copy the
  reference's grid + converged potential, swap in the new detector + impurity
  model, and re-converge via `update_till_convergence!` (no refinement, no
  cold restart). Much faster when `det` is a small perturbation of the
  reference. Then applies
  [`adapt_potential_to_depletion_plus_offset!`](@ref) so `det.Vop` lands at
  `Vdep + 500 V` rather than inheriting the reference's bias — the
  superposition step that dominates compute time here, but kept for
  consistent `Vop` semantics with the other `characterize!` paths.
"""
characterize!(det::DetectorDesign{T}, boule::CrystallineBoule{T}, args...; kwargs...) where {T<:AbstractFloat} =
    characterize!(det, impurity_model_from_boule(boule, det), args...; kwargs...)

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
        end
        if length(refinement_limits) <= 2 || is_depleted(sim.point_types)
            adapt_potential_to_depletion_plus_offset!(det, sim; verbose = verbose)
            populate_design!(det, sim, det.Vop; verbose = verbose)
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
    verbose::Bool = false
) where {T<:AbstractFloat}
    sim = deepcopy(reference_simulation) #should we edit simulation in place or deepcopy?
    sim.detector = SolidStateDetector{T}(det, imp_model)
    
    update_till_convergence!( sim, ElectricPotential,
                                    n_iterations_between_checks = 1000,
                                    max_n_iterations = 50000,
                                    depletion_handling = true,
                                    verbose = verbose
                            )
    mark_bulk_bits!(sim.point_types.data)
    mark_undep_bits!(sim.point_types.data, sim.imp_scale.data)
    det.is_simulated = true

    initialize_design!(det)
    if is_depleted(sim.point_types)
        # TODO: if perturbation is small enough, is this necessary? It's adding most of the
        # compute time. Skipping it would just compute a new depletion voltage but Vop would
        # be inherited from the reference instead of landing at exactly Vdep + offset.
        adapt_potential_to_depletion_plus_offset!(det, sim; verbose = verbose)
        populate_design!(det, sim, det.Vop; verbose = verbose)
    end
    sim
end
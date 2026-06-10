# Tutorial

This walkthrough builds and Characterizes a generic inverted-coaxial point-contact
(ICPC) germanium detector cut from a synthetic boule. All numbers below are made
up — the goal is to show the *workflow*, not reproduce a specific detector.

The full pipeline:

1. Define a [`CrystallineBoule`](@ref) from impurity measurements + an external profile.
2. Fit the impurity-density curve.
3. Define an [`InvertedCoaxDesign`](@ref).
4. Characterize it with [`characterize!`](@ref) — runs the SSD field solver and
   populates depletion voltage, minimum bulk E-field, mass.
5. Inspect and plot.
6. Iterate on the geometry without re-doing the cold solve.

## Setup

```julia
using LegendDetectorDesign
const LDD = LegendDetectorDesign

using Unitful
using SolidStateDetectors
using Plots
using LsqFit
```

`T` controls the floating-point precision used throughout the geometry, the
impurity model parameters, and the SSD grid. `Float32` is the standard choice
and matches LEGEND production simulations:

```julia
T = Float32
```

## 1. Define a crystalline boule

A [`CrystallineBoule`](@ref) bundles the boule's identity, its external profile
(a [`BouleGeometry`](@ref)), the raw impurity measurements, and the parameters
of an impurity-density fit. Start with the measured Hall samples and a
piecewise radius profile:

```julia
z_hall        = [ 0,  30,  60,  90, 120] * u"mm"
impurity_hall = [0.4, 0.5, 0.7, 1.0, 1.5] * u"1e10/cm^3"

# axially-symmetric boule outline; measured radius at three z stations
z_radius     = [0, 65, 130] * u"mm"
r_radius     = [42, 42, 41] * u"mm"

boule = CrystallineBoule(T;
    name             = "EXBoule",
    order            = "01",
    impurity_model   = :linear_exponential_boule,
    impurity_hall    = impurity_hall,
    z_hall           = z_hall,
    geometry         = BouleGeometry(T; z = z_radius, radius = r_radius),
)
```

[`BouleGeometry`](@ref) accepts non-uniform `z` samples and fits a
Fritsch–Carlson monotonic cubic spline through them — the same interpolator
SSD uses for `SplineBouleImpurityDensity`.

The four built-in impurity-model symbols accepted by `impurity_model =` are:

| symbol                          | formula                                         |
|---------------------------------|-------------------------------------------------|
| `:linear_boule`                 | `a + b·z`                                       |
| `:parabolic_boule`              | `a + b·z + c·z²`                                |
| `:linear_exponential_boule`     | `a + b·z + n·exp((z − l)/m)`                    |
| `:parabolic_exponential_boule`  | `a + b·z + c·z² + n·exp((z − l)/m)`             |

The exponential-tail families are useful when the resistivity rises sharply
near one end of the boule.

## 2. Fit the impurity model

Use any nonlinear-least-squares tool — LsqFit.jl works directly with the curve
functions returned by [`fit_function`](@ref):

```julia
fit = curve_fit(
    boule.impurity_model,        # = LDD.linear_exponential_boule
    boule.z_hall,
    boule.impurity_hall,
    [1.0, 0.1, 1.0, 70.0, 30.0], # initial guess for [a, b, n, l, m]
)

boule.impurity_model_parameters = fit.param
```

`boule.impurity_model_parameters` is stored as bare numbers in internal units.
To recover the unitful form, use [`get_unitful_property`](@ref):

```julia
get_unitful_property(boule, :impurity_model_parameters)
```

## 3. Visualise the impurity profile

[`get_impurity_density`](@ref) evaluates the fitted curve at any axial
position, broadcasting over arrays:

```julia
x_dense = range(-5u"mm", maximum(z_hall) + 5u"mm", length = 200)

scatter(z_hall, impurity_hall,
    label   = "Hall samples",
    xguide  = "z along boule",
    yguide  = "Impurity density",
    frame   = :box,
)
plot!(x_dense, get_impurity_density(boule, x_dense),
    label = "fit", lw = 2,
)
```

## 4. Define the detector

[`InvertedCoaxDesign`](@ref) builds an ICPC geometry from named dimensions and
wraps it in a [`DetectorDesign`](@ref). All length kwargs accept `Number`
(interpreted as mm) or a `Unitful.Quantity`; angles default to degrees.

```julia
det = InvertedCoaxDesign(T;
    name                   = "ExampleICPC",
    height                 = 80,
    radius                 = 40,
    pc_radius              = 7,
    borehole_pc_gap        = 25,
    borehole_radius        = 5,
    borehole_taper_height  = 35,
    borehole_taper_angle   = 5,
    top_taper_height       = 1,  
    top_taper_angle        = 45,
    bottom_taper_height    = 1,
    bottom_taper_angle     = 45,
    dead_layer_depth       = 0.8,
    offset                 = 105,   # where the p⁺ contact sits on the boule axis
)
```

The constructor tags the geometry [`ValidGeometry`](@ref) only when the
dimensions satisfy basic consistency checks *and* `offset ≥ height` (so the
detector fits inside the boule). Invalid designs return with `mass = missing`
and won't simulate.

`det.offset` is the position of the p⁺ point contact relative to the boule's
seed end — it tells the impurity model which slice of the boule this crystal
was cut from.

## 5. Plot the detector

The plot extension loads automatically when `Plots` is present. For a quick
look at the geometry alone:

```julia
plot(det.geometry)                                # filled cross-section
plot(det.geometry, include_measurements = true)   # with dimensions called out
```

For a full measured drawing suitable for cross-checking against a
manufacturer's spec sheet — calls out every dimension, the p⁺ spot pattern,
and embeds title / order labels:

```julia
plot(det, technical_drawing = true, crystal_prefix = "EX", corner_rounding = :both)
```

`technical_drawing = true` uses the longer-form recipe for `DetectorDesign`;
`include_measurements = true` is the lighter equivalent for the bare geometry.

## 6. Plot the detector in the boule

The package ships a two-argument recipe `plot(boule, det)` that draws the
detector's silhouette positioned at `det.offset` along the boule profile.
This is the cut-position visualisation you want when deciding where to slice
a new design from a known boule:

```julia
plot(boule, det)                            # boule outline + detector at offset
plot(boule, det, technical_drawing = true)  # full layout with section views A-A / B-B
```

In the technical-drawing form the recipe also adds circular section cuts
through the boule at the detector's top and bottom planes plus offset
markers — useful for verifying the boule has enough radius to accommodate
the detector everywhere along its height.

## 7. Characterize the detector

[`characterize!`](@ref) runs the SSD field solver and fills in `det.Vdep`,
`det.Vop`, `det.Emin`, `det.Emin_pos`, `det.mass`:

```julia
sim = characterize!(det, boule,
    refinement_limits = [0.2, 0.1, 0.05, 0.02],
    verbose           = true,
)
det     # prints a unicode summary
```

What happens internally:

1. The impurity model is built from `boule.impurity_model_parameters`,
   offset-adjusted by `det.offset`.
2. SSD is run at the maximum allowed voltage (`Vmax`, default `5000 V`) with
   two coarse refinement passes.
3. If the detector depletes, the grid is refined further and the converged
   potential is projected to `Vdep + 500 V` analytically using the
   weighting-potential superposition trick — no second cold solve required.
4. `populate_design!` extracts the depletion voltage, the minimum bulk
   E-field, and its `(r, z)` location.

The returned `sim::SSD.Simulation` is kept around so you can plot or
post-process it.

## 6. Inspect the results

Each detector field is a bare number in internal units:

```julia
det.Vdep              # depletion voltage (V)
det.Vop               # operating voltage used (V)
det.Emin              # minimum bulk E field (V/cm)
det.Emin_pos          # (r, z) of that minimum (mm, mm)
det.mass              # crystal mass (g)
```

Convert to unitful values when you need them — e.g.
`det.Vdep * u"V"`. Convenience accessors:

```julia
get_unitful_property(det, :offset)   # det.offset in mm
```

## 7. Plot the detector

The plot extension is loaded automatically when `Plots` is present:

```julia
plot(det, technical_drawing = true)
```

`technical_drawing = true` produces a measured drawing with all the
dimensions called out — useful for ordering / cross-checking against a
manufacturer's spec sheet. For a lighter-weight schematic, use
`include_measurements = true` on the geometry directly:

```julia
plot(det.geometry, include_measurements = true)
```

## 8. Plot the electric field

`sim.electric_field` and `sim.electric_potential` are SSD scalar/vector fields
on the SSD grid — they plot directly:

```julia
plot(sim.electric_field,
    full_det               = true,
    color                  = cgrad(:inferno, scale = :exp),
    clims                  = (0, 3) .* u"kV/cm",
    contours_equal_potential = true,
    levels                 = 15,
    linecolor              = :white,
)
```

For just the depletion/point-type map:

```julia
plot(sim.point_types, full_det = true)
```

## 9. Iterate on the geometry — warm-start

You'll often want to sweep one dimension while keeping the rest fixed. The
keyword-only clone-constructor [`InvertedCoaxGeometry`](@ref)`(geo; …)` does
exactly that:

```julia
geo2 = InvertedCoaxGeometry(det.geometry;
    pc_radius = T(5.5),
    dead_layer_depth = T(0.5),
)
det2 = InvertedCoaxDesign(det, geo2)
```

The previous simulation is a good initial guess for the new one — use the
**warm-start** form of `characterize!` that takes the prior `Simulation` and
re-converges from its potential instead of cold-resetting:

```julia
characterize!(det2, boule, sim)
det2
```

SSD continues from the existing grid and potential via
`update_till_convergence!`, typically converging much faster than a fresh
cold solve. Equally valid when sweeping the bias voltage on an unchanged
geometry — just pass the new `Vop`.

## 10. Persist metadata

Three companion functions serialise to LEGEND-format metadata containers:

```julia
geo_meta    = LDD.geo_to_meta(det.geometry; Vop = det.Vop, name = det.name)
det_meta    = LDD.design_to_meta(det)        # PropDict, suitable for SSD round-trip
boule_meta  = LDD.boule_to_meta(boule, det)  # OrderedDict, LEGEND boule YAML layout
```

None of them write to disk — choose the writer yourself:

```julia
using YAML
YAML.write_file("boule_$(boule.order).yaml", boule_meta)
```

The inverse `meta_to_geo` / `meta_to_design` parse a `PropDict` back into a
geometry / design.

## Next steps

- See the [API reference](api.md) for the full list of types, exported
  functions, and the four available impurity models.
- For lower-level SSD operations — direct potential calculation, charge drift
  simulation, weighting-potential extraction — consult the
  [SolidStateDetectors.jl](https://github.com/JuliaPhysics/SolidStateDetectors.jl)
  documentation; everything in this package returns standard SSD objects.

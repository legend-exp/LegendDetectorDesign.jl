# Tutorial

This walkthrough builds and characterizes a generic inverted-coaxial point-contact
(ICPC) germanium detector cut from a boule. All numbers below are
placeholders — the goal is to show the *workflow*, not reproduce a specific detector.

The full pipeline:

1. Define a [`CrystallineBoule`](@ref) from impurity measurements + an external profile.
2. Fit the impurity-density curve.
3. Define an [`InvertedCoaxDesign`](@ref).
4. Plot the detector and its position in the boule.
5. Characterize it with [`characterize!`](@ref) — runs the SSD field solver and
   populates depletion voltage, minimum bulk E-field, mass.
6. Inspect and plot the field.
7. Iterate on the geometry without re-doing the cold solve.
8. Serialize to LEGEND-format metadata containers.

## Setup

```@example tutorial
ENV["GKSwstype"] = "100" # hide

using LegendDetectorDesign
const LDD = LegendDetectorDesign

using Unitful
using SolidStateDetectors
using Plots
using LsqFit
```

`T` controls the floating-point precision used throughout the geometry, the
impurity-model parameters, and the SSD grid. `Float32` is the standard choice
and matches LEGEND production simulations:

```@example tutorial
T = Float32
```

## 1. Define a crystalline boule

A [`CrystallineBoule`](@ref) bundles the boule's identity, its external profile
(a [`BouleGeometry`](@ref)), the raw impurity measurements, and the parameters
of an impurity-density fit. Start with the measured Hall samples and a
piecewise radius profile:

```@example tutorial
z_hall        = [ 0, 25, 55, 85, 110] * u"mm"
impurity_hall = [0.40, 0.55, 0.70, 0.85, 1.05] * u"1e10/cm^3"

# axially-symmetric boule outline; measured radius at three z stations
z_radius = [0, 60, 120] * u"mm"
r_radius = [42, 42, 41] * u"mm"

boule = CrystallineBoule(T;
    name           = "EXBoule",
    order          = "01",
    impurity_model = :linear_exponential_boule,
    impurity_hall  = impurity_hall,
    z_hall         = z_hall,
    geometry       = BouleGeometry(T; z = z_radius, radius = r_radius),
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

```@example tutorial
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

```@example tutorial
get_unitful_property(boule, :impurity_model_parameters)
```

## 3. Visualize the impurity profile

[`get_impurity_density`](@ref) evaluates the fitted curve at any axial
position, broadcasting over arrays:

```@example tutorial
x_dense = range(-5u"mm", maximum(z_hall) + 5u"mm", length = 200)

scatter(z_hall, impurity_hall;
    label  = "Hall samples",
    xguide = "z along boule",
    yguide = "Impurity density",
    frame  = :box,
    size   = (640, 360),
)
plot!(x_dense, get_impurity_density(boule, x_dense), label = "fit", lw = 2)
```

## 4. Define the detector

[`InvertedCoaxDesign`](@ref) builds an ICPC geometry from named dimensions and
wraps it in a [`DetectorDesign`](@ref). All length kwargs accept `Number`
(interpreted as mm) or a `Unitful.Quantity`; angles default to degrees.

```@example tutorial
det = InvertedCoaxDesign(T;
    name                   = "ExampleICPC",
    height                 = 80,
    radius                 = 38,
    pc_radius              = 8,
    borehole_pc_gap        = 30,
    borehole_radius        = 5,
    borehole_taper_height  = 40,
    borehole_taper_angle   = 5,
    top_taper_height       = 1,
    top_taper_angle        = 45,
    bottom_taper_height    = 1,
    bottom_taper_angle     = 45,
    dead_layer_depth       = 0.8,
    offset                 = 100,
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

```@example tutorial
plot(det.geometry, size = (300, 500))
```

The same recipe with `include_measurements = true` calls out every dimension:

```@example tutorial
plot(det.geometry, include_measurements = true, size = (500, 700))
```

For a full measured drawing — the kind you'd attach to an order to a
manufacturer — the `DetectorDesign` recipe adds the p⁺ spot pattern, title,
and order labels:

```@example tutorial
plot(det, technical_drawing = true, crystal_prefix = "EX", corner_rounding = :both)
```

## 6. Plot the detector in the boule

The package ships a two-argument recipe `plot(boule, det)` that draws the
detector's silhouette positioned at `det.offset` along the boule profile —
the cut-position visualization you want when picking where to slice a new
design from a known boule:

```@example tutorial
plot(boule, det, size = (1100, 360))
```

`technical_drawing = true` adds circular section cuts through the boule at
the detector's top and bottom planes, plus offset markers — useful for
verifying the boule has enough radius to accommodate the detector everywhere
along its height:

```@example tutorial
plot(boule, det, technical_drawing = true, size = (1100, 400))
```

## 7. Characterize the detector

[`characterize!`](@ref) runs the SSD field solver and fills in `det.Vdep`,
`det.Vop`, `det.Emin`, `det.Emin_pos`, `det.mass`. Heavier refinement gives
finer field maps at the cost of build time; the tutorial uses a moderate
two-pass schedule:

```@example tutorial
sim = characterize!(det, boule, refinement_limits = [0.2, 0.1, 0.05])
det
```

What happens internally:

1. The impurity model is built from `boule.impurity_model_parameters`,
   offset-adjusted by `det.offset`.
2. SSD is run at the maximum allowed voltage (`Vmax`, default `5000 V`) with
   two coarse refinement passes.
3. If the detector depletes, the converged potential is projected to
   `Vdep + 500 V` analytically using the weighting-potential superposition
   trick — no second cold solve required.
4. `populate_design!` extracts the depletion voltage, the minimum bulk
   E-field, and its `(r, z)` location.

The returned `sim::SSD.Simulation` is kept around so you can plot or
post-process it.

## 8. Inspect the results

Each detector field is a bare number in internal units:

```@example tutorial
(Vdep = det.Vdep, Vop = det.Vop, Emin = det.Emin, Emin_pos = det.Emin_pos, mass = det.mass)
```

Convert to unitful values when you need them — e.g. `det.Vdep * u"V"`.
Convenience accessor for the cutting offset:

```@example tutorial
get_unitful_property(det, :offset)
```

## 9. Plot the electric field

`sim.electric_field` and `sim.electric_potential` are SSD scalar / vector
fields on the SSD grid — they plot directly:

```@example tutorial
plot(sim.electric_field;
    full_det       = true,
    color          = cgrad(:inferno, scale = :exp),
    xunit          = u"mm",
    yunit          = u"mm",
    zunit          = u"kV/cm",
    clims          = (0, 3),
    legendfontsize = 6,
)
plot_electric_fieldlines!(sim; sampling = 5u"mm", lw = 0.5, title = "")
plot!(sim.detector.contacts[1]; st = :slice, lw = 1)
plot!(sim.detector.contacts[2]; st = :slice, c = :lightgrey, lw = 2)
```

The depletion / point-type map is one plot call:

```@example tutorial
plot(sim.point_types, full_det = true)
```

## 10. Iterate on the geometry — warm-start

You'll often want to sweep one dimension while keeping the rest fixed. The
keyword-only clone-constructor [`InvertedCoaxGeometry`](@ref)`(geo; …)` does
exactly that:

```@example tutorial
geo2 = InvertedCoaxGeometry(det.geometry;
    pc_radius        = T(8.5)
)
det2 = InvertedCoaxDesign(det, geo2)
```

The previous simulation is a good initial guess for the new one — use the
**warm-start** form of `characterize!` that takes the prior `Simulation` and
re-converges from its potential instead of cold-resetting:

```@example tutorial
characterize!(det2, boule, sim)
det2
```

SSD continues from the existing grid and potential via
`update_till_convergence!`, typically converging much faster than a fresh
cold solve.

## 11. Persist metadata

Two companion functions serialize to LEGEND-format metadata containers:

```@example tutorial
LDD.design_to_meta(det)
```

```@example tutorial
LDD.boule_to_meta(boule, det)
```

`design_to_meta` is the `DetectorDesign` flavor (it just delegates to `geo_to_meta` with the design's `Vop` and `name`). None of them write
to disk — choose the writer yourself:

```julia
using YAML
YAML.write_file("boule_$(boule.order).yaml", LDD.boule_to_meta(boule, det))
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

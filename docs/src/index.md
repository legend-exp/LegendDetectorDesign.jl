# LegendDetectorDesign.jl

Tools for designing and simulating high-purity germanium (HPGe)
detectors for the [LEGEND experiment](https://legend-exp.org/). Built on
top of [`SolidStateDetectors.jl`](https://github.com/JuliaPhysics/SolidStateDetectors.jl)
and [`LegendDataManagement.jl`](https://github.com/legend-exp/LegendDataManagement.jl).

## What it does

The package lets you go from raw boule measurements to a simulated detector
with a few high-level types:

| object | role |
|---|---|
| [`CrystallineBoule`](@ref) | A grown germanium boule with its external profile, fit-ready impurity model, and raw Hall / resistivity measurements. |
| [`BouleGeometry`](@ref) | Axially-symmetric boule outline; monotonic spline through measured radii (matches SSD's `SplineBouleImpurityDensity` idiom). |
| [`InvertedCoaxGeometry`](@ref) | ICPC detector geometry parametrized by the usual `(height, radius, pc_radius, borehole_*, *_taper_*, dead_layer_depth)` set. |
| [`DetectorDesign`](@ref) | A geometry plus the simulated characterization results (depletion voltage, operating voltage, minimum bulk field, mass). |
| [`characterize!`](@ref) | Runs the SSD field solver on a design + impurity model and populates the depletion/operating/field/mass slots. |

Plot recipes for every type are loaded automatically when `Plots` is present,
so cross-sections, technical drawings, and field maps are a single `plot(...)`
call.

## Installation

`LegendDetectorDesign.jl` lives in the
[`legend-julia-registry`](https://github.com/legend-exp/legend-julia-registry)
along with the rest of the LEGEND Julia stack. Add the registry once, then
install the package the usual way:

```julia
using Pkg
pkg"registry add https://github.com/JuliaRegistries/General"
pkg"registry add https://github.com/legend-exp/legend-julia-registry"
pkg"add LegendDetectorDesign"
```

## Quick start

A complete walkthrough — boule → fit → detector → field solver → plots — is
in the [tutorial](tutorial.md). A short teaser:

```julia
using LegendDetectorDesign, Unitful, Plots

T = Float32

boule = CrystallineBoule(T;
    name = "EXBoule", order = "01",
    impurity_model = :linear_exponential_boule,
    impurity_hall  = [0.40, 0.55, 0.70, 0.85, 1.05] * u"1e10/cm^3",
    z_hall         = [   0,   25,   55,   85,  110] * u"mm",
    geometry       = BouleGeometry(T; z = [0, 60, 120] * u"mm", radius = [42, 42, 41] * u"mm"),
)

det = InvertedCoaxDesign(T;
    name = "ExampleICPC",
    height = 80, radius = 38, pc_radius = 8,
    borehole_pc_gap = 30, borehole_radius = 5,
    borehole_taper_height = 40, dead_layer_depth = 0.8, offset = 100,
)

# Fit the impurity model, run the field solver, plot the detector + field
# (see the tutorial for the full example)
```

## Where to read next

- **[Tutorial](tutorial.md)** — end-to-end synthetic example with rendered
  plots and field maps.
- **[API reference](api.md)** — every exported type, constructor, and helper.
- **[LICENSE](LICENSE.md)** — MIT.

## Related packages

- [`SolidStateDetectors.jl`](https://juliaphysics.github.io/SolidStateDetectors.jl/stable/) —
  the underlying field solver. `LegendDetectorDesign` returns standard SSD
  objects (`Simulation`, `SolidStateDetector`, `ElectricField`, …); anything
  you can do directly in SSD works on the results from `characterize!`.
- [`LegendDataManagement.jl`](https://github.com/legend-exp/LegendDataManagement.jl) —
  reads / writes the LEGEND metadata formats this package targets via
  `geo_to_meta` / `design_to_meta` / `boule_to_meta`.

## Acknowledgements

Developed for the LEGEND collaboration by David Hervas Aguilar. Contributions, issues, and design
discussions are welcome on
[GitHub](https://github.com/legend-exp/LegendDetectorDesign.jl).

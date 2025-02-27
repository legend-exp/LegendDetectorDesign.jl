# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using LegendDetectorDesign

# Doctest setup
DocMeta.setdocmeta!(
    LegendDetectorDesign,
    :DocTestSetup,
    :(using LegendDetectorDesign);
    recursive=true,
)

makedocs(
    sitename = "LegendDetectorDesign",
    modules = [LegendDetectorDesign],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://legend-exp.github.io/LegendDetectorDesign.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    warnonly = ("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/legend-exp/LegendDetectorDesign.jl.git",
    forcepush = true,
    push_preview = true,
)

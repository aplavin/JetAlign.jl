module JetAlign

using DataManipulation
using StaticArrays
using AccessorsExtra
using InterferometricModels
using InterferometricModels: fwhm_to_σ
using VLBIFiles
using AxisKeysExtra
using AxisKeysExtra.AxisKeys: keyless_unname
using StatsBase: geomean
using RectiGrids
using IntervalSets
using Unitful, UnitfulAngles
using Statistics
using LinearAlgebra: norm
using PyFormattedStrings
using Colors
import ColorSchemes
using SQLite, SQLCollections
using Scratch
using UUIDs
import TOML

import CImGui as ig
import ImPlot
import ImGuiThemes
import GLFW, ModernGL   # activate CImGui's GlfwOpenGL3 backend
using ImPlotExtra: image!, contour!, colorbar!

export imalign_interactive

include("utils.jl")
include("string_alignment.jl")
include("image_processing.jl")
include("clean_components.jl")
include("app.jl")
include("ui.jl")


"""
    imalign_interactive(files; shiftsfile=nothing, colors=(colorant"red", colorant"blue"), persiststate=true)

Open the interactive alignment window for two images (FITS files with clean components) and return when it
is closed. The lower-frequency image is the one shifted.

- `shiftsfile`: a shifts TOML written by "Save shifts" for the same `files`; its shifts become the multishift points.
- `colors`: the two images' colors, by increasing frequency.
- `persiststate`: restore the UI state of this file pair on launch and save it while running.
"""
imalign_interactive(files; kwargs...) = run!(App(files; kwargs...))


using PrecompileTools

@compile_workload let
    map_fits = joinpath(pkgdir(VLBIFiles), "test", "data", "map.fits")
    app = App([map_fits, map_fits]; persiststate=false)
    compute_maps(app.images, maps_inputs(app))
    set_shifts!(app, [SVector(0.1, 0.2), SVector(-0.3, 0.1)])
    maps = compute_maps(app.images, maps_inputs(app))
    map(int -> si_color(app.ui.sicolorrange, int), maps.si)
    compute_xcorr(app.images, app.xcorr_halfwidth, xcorr_inputs(app))
end

end

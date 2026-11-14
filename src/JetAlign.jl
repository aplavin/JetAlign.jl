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
using ArgParse

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


"Command-line interface of the `jetalign` app."
function cli_settings()
    s = ArgParseSettings(prog="jetalign", description="Interactively align two multi-frequency VLBI images.")
    @add_arg_table! s begin
        "--shifts"
            help = "multishift points from a file written by \"Save shifts\" for the same files"
            metavar = "FILE"
        "--no-persist"
            help = "don't restore or save the UI state for this file pair"
            action = :store_true
        "image1"
            help = "FITS image with clean components"
            required = true
        "image2"
            help = "FITS image of the same source at another frequency"
            required = true
    end
    s
end

function (@main)(argv)
    args = parse_args(argv, cli_settings())
    imalign_interactive([args["image1"], args["image2"]]; shiftsfile=args["shifts"], persiststate=!args["no-persist"])
    return 0
end


using PrecompileTools

@compile_workload let
    parse_args(["--shifts", "s.toml", "a.fits", "b.fits"], cli_settings())
    map_fits = joinpath(pkgdir(VLBIFiles), "test", "data", "map.fits")
    app = App([map_fits, map_fits]; persiststate=false)
    compute_maps(app.images, maps_inputs(app))
    set_shifts!(app, [SVector(0.1, 0.2), SVector(-0.3, 0.1)])
    maps = compute_maps(app.images, maps_inputs(app))
    map(int -> si_color(app.ui.sicolorrange, int), maps.si)
    compute_xcorr(app.images, app.xcorr_halfwidth, xcorr_inputs(app))
end

end

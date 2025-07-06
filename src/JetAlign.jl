module JetAlign

using DataManipulation
using StaticArrays
using AccessorsExtra
using InterferometricModels
using VLBIFiles; using VLBIFiles: coords, flux
using AxisKeysExtra
using StatsBase
using RectiGrids
using VLBIPlots
using IntervalSets
using Unitful, UnitfulAngles, UnitfulAstro
using Statistics
using LinearAlgebra: norm
using ImageFiltering
using PyFormattedStrings
using MakieExtra
import GLMakie
using SQLite, SQLCollections
using Scratch

export imalign_interactive

include("utils.jl")
include("string_alignment.jl")
include("observables.jl")
include("image_processing.jl")
include("clean_components.jl")

const RAN_ONCE = Ref(false)

function imalign_interactive(files; colors=[:red, :blue], display=true)
    length(files) == 2 || error("Exactly two files must be provided for alignment, got $(length(files))")
    files = map(expanduser, files)
    if !all(isfile, files)
        error("Files not found:\n$(@p files filter(!isfile) join(__, '\n'))")
    end
    RAN_ONCE[] || (println("running JetAlign for the first time, this may take a while..."); RAN_ONCE[] = true)

    images_with_data = @p map(files) do f
        fimg = VLBI.load(f, read_data=false)
        mod = VLBI.load(MultiComponentModel, f)
        mod = @modify(ustrip, mod |> RecursiveOfType(ModelComponent) |> flux)
        (;freq=VLBI.frequency(fimg), mod, beam=VLBI.load(Beam, f))
    end sort(by=_.freq)
    avgcircbeam = @p images_with_data map(fwhm_average(_.beam)) geomean() construct(Beam{CircularGaussian}, fwhm_average => __ |> u"mas")

    radecs = @p range(0u"mas" ± 10*fwhm_average(avgcircbeam), length=512) grid(ra=reverse(__), dec=__)

    fig = Figure(size=(1300, 500))

    σobs = Observable(ustrip(u"mas", avgcircbeam.σ))
    curbeam = @lift beam(CircularGaussian, σ=$σobs)

    shiftobs = Observable(SVector(0., 0.))
    img_a = @lift cc_to_image(images_with_data[1].mod; σ=$σobs*u"mas", radecs, shift=SVector($shiftobs)u"mas")
    img_b = @lift cc_to_image(images_with_data[2].mod; σ=$σobs*u"mas", radecs)
    SI = @lift spectral_index($img_a, $img_b, map(:freq, images_with_data); threshold=0.3e-3)

    shiftmax = 1.5*fwhm_average(avgcircbeam) |> u"mas"

    corrmasksize = Observable(1.)
    xcorr_unsm = @lift unsmoothed_xcorr_ccs(CorrelationKind.Raw(),
        mask_ccs(images_with_data[2].mod, $corrmasksize*u"mas"),
        mask_ccs(images_with_data[1].mod, $corrmasksize*u"mas");
        radecs=(@p range(0u"mas"±2shiftmax, length=256) grid(ra=reverse(__), dec=__)));
    xcorrimg = @lift smooth($xcorr_unsm; σ=$σobs*u"mas")

    curpos = Observable(SVector(NaN, NaN))
    curpos_fplt = @lift FPlot($curpos, identity; color=(:white, 0.5), strokewidth=1, strokecolor=Ref(:black), marker=marker_lw(:cross, 0.5), markersize=15)

    shiftax = Axis(fig[1,1], limits=(-1, 1, -1, 1) .* ustrip(u"mas", shiftmax), aspect=DataAspect(), xreversed=true, xlabel="ra (mas)", ylabel="dec (mas)", backgroundcolor=:black)
    image!(xcorrimg, colormap=:turbo, colorrange=@lift (0.6, 1) .* maximum(abs, $xcorrimg))
    text!((0.99, 0.99), text="Color = X-corr\nImage shift:\ndrag with the mouse", space=:relative, align=(:right, :top), color=:white)
    vlines!(0, color=:white)
    hlines!(0, color=:white)
    scatter!(shiftobs; color=:white, strokecolor=:black, strokewidth=1)
    beampoly!(curbeam, color=(:white, 0), strokecolor=(:white, 0.3), strokewidth=1)

    ax_si, plt_si = image(fig[1,3], SI, colormap=:turbo, colorrange=(-1.8, 1.2), axis=(xlabel="", ylabel=""), interpolate=true)
    text!((0.99, 0.99), text="Spectral index", space=:relative, align=(:right, :top))
    scatter!(curpos_fplt)
    Colorbar(fig[1,4], plt_si; height=250, tellheight=false)

    ax, _ = contour(fig[1, 2], img_a, levels=(@p maximum(img_a[]) maprange(log, __/1e3, __*0.9, length=15)), color=colors[1])
    contour!(img_b, levels=(@p maximum(img_b[]) maprange(log, __/1e3, __*0.9, length=15)), color=colors[2])
    poly!((@lift Circle(Makie.Point(0.,0), $corrmasksize)), color=(:black, 0.12), strokewidth=0.4, strokecolor=:black)
    text!((0.01, 0.99), text="Stokes I", space=:relative, align=(:left, :top))
    text!((0.99, 0.99), text=(
        rich(f"{images_with_data[1].freq:.1f}: {maximum(img_a[]):.2f} Jy/bm", color=colors[1])
        *"\n"*
        rich(f"{images_with_data[2].freq:.1f}: {maximum(img_b[]):.2f} Jy/bm", color=colors[2])
    ), space=:relative, align=(:right, :top))
    beampoly!(ax, curbeam, color=(:black, 0), strokecolor=:black, strokewidth=1)
    beampoly!(images_with_data[1].beam|>ustrip, color=(:white, 0), strokecolor=(colors[1], 0.3), strokewidth=1)
    beampoly!(images_with_data[2].beam|>ustrip, color=(:white, 0), strokecolor=(colors[2], 0.3), strokewidth=1)
    scatter!(curpos_fplt)

    linkaxes!(ax_si, ax)


    map!(identity, shiftobs, mouse_position_obs(shiftax; key=Exclusively(Mouse.left)); update=false)
    map!(curpos, mouse_position_obs(ax; hold=false), mouse_position_obs(ax_si; hold=false)) do p1, p2
        all(!isnan, p1) ? p1 : p2
    end

    sg = SliderGrid(fig[0,:][1,2],
        (label = "Beam FWHM", range=0.05:0.01:1.5, format = "{:.2f} mas", startvalue=InterferometricModels.σ_to_fwhm(σobs[])),
        (label = "X-corr mask radius", range=0:0.1:3, format = "{:.2f} mas", startvalue=corrmasksize[]),
        tellheight=true, tellwidth=false)
    map!(InterferometricModels.fwhm_to_σ, σobs, sg.sliders[1].value)
    map!(identity, corrmasksize, sg.sliders[2].value; update=true)

    let hl_strs = highlight_string_differences(files; colors)
        Label(fig[0,:][1,1], rich(hl_strs[1] * "\n" * hl_strs[2], color=:gray70, font="JuliaMono"))
    end
    colorbuffer(fig)  # to ensure default finallimits are set

    statekey = @NamedTuple{f1::String, f2::String}(NTuple{2}(files))
    stateobs = (;
        slider_ixs = map(sl -> sl.selected_index, sg.sliders),
        shiftobs,
        ax_limits = obs_getset(ax.finallimits, ax.limits),
        shiftax_limits = obs_getset(shiftax.finallimits, shiftax.limits),
    )
    stateobs_single = lift_getset(
        () -> modify(o -> o[], stateobs, RecursiveOfType(Observable)),
        state -> modify((o, v) -> o[] = v, stateobs, RecursiveOfType(Observable), state),
        getall(stateobs, RecursiveOfType(Observable))...
    )
    statedefault = stateobs_single[]
    on(Button(fig[0,:][1,3], label="Reset").clicks) do _
        stateobs_single[] = statedefault
    end
    setup_state_persistence(statekey, stateobs_single)

    if display
        MakieExtra.show_gl_icon_in_dock()
        Base.display(fig)
    else
        return fig
    end
end

end

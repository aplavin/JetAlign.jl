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
include("image_processing.jl")
include("correlation.jl")

function imalign_interactive(files; colors=[:red, :blue], display=true)
    length(files) == 2 || error("Exactly two files must be provided for alignment, got $(length(files))")
    files = map(expanduser, files)
    if !all(isfile, files)
        error("Files not found:\n$(@p files filter(!isfile) join(__, '\n'))")
    end

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
        radecs=(@p range(0u"mas"±shiftmax, length=128) grid(ra=reverse(__), dec=__)));
    xcorrimg = @lift smooth($xcorr_unsm; σ=$σobs*u"mas")

    curpos = Observable(SVector(NaN, NaN))
    curpos_fplt = @lift FPlot($curpos, identity; color=(:white, 0.5), strokewidth=1, strokecolor=Ref(:black), marker=marker_lw(:cross, 0.5), markersize=15)

    shiftax = Axis(fig[1,1], limits=(-1, 1, -1, 1) .* ustrip(u"mas", shiftmax), aspect=DataAspect(), xreversed=true, xlabel="ra (mas)", ylabel="dec (mas)", backgroundcolor=:black)
    image!(xcorrimg, colormap=:turbo, colorrange=@lift (0.6, 1) .* maximum(abs, $xcorrimg))
    text!((0.99, 0.99), text="Image shift\nColor = X-corr", space=:relative, align=(:right, :top), color=:white)
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


    on(events(shiftax).mousebutton, priority=10) do evt
        if is_mouseinside(shiftax) && evt.button == Makie.Mouse.left
            shiftobs[] = Tuple(mouseposition(shiftax))
            Consume()
        end
    end

    on(events(fig).mouseposition, priority=10) do pos
        if is_mouseinside(shiftax) && ispressed(fig, Makie.Mouse.left)
            shiftobs[] = Tuple(mouseposition(shiftax))
            return Consume()
        end
    end

    on(events(fig).mouseposition) do pos
        if is_mouseinside(ax) || is_mouseinside(ax_si)
            curpos[] = Tuple(mouseposition(is_mouseinside(ax) ? ax : ax_si))
        else
            curpos[] = (NaN, NaN)
        end
    end

    sg = SliderGrid(fig[0,:][1,2],
        (label = "Beam FWHM", range=0.05:0.01:1.5, format = "{:.2f} mas", startvalue=InterferometricModels.σ_to_fwhm(σobs[])),
        (label = "X-corr mask radius", range=0:0.1:3, format = "{:.2f} mas", startvalue=corrmasksize[]),
        tellheight=true, tellwidth=false)
    on(sg.sliders[1].value) do fwhm
        σobs[] = InterferometricModels.fwhm_to_σ(fwhm)
    end
    notify(sg.sliders[1].value)
    map!(identity, corrmasksize, sg.sliders[2].value)
    notify(sg.sliders[2].value)

    let hl_strs = highlight_string_differences(files; colors)
        Label(fig[0,:][1,1], rich(hl_strs[1] * "\n" * hl_strs[2], color=:gray70, font="JuliaMono"))
    end
    colorbuffer(fig)  # to ensure default finallimits are set

    statekey = NamedTuple{(:f1, :f2)}(NTuple{2}(files))
    stateobs = (;
        slider_ixs = map(sl -> sl.selected_index, sg.sliders),
        shiftobs,
        ax_limits = obs_getset(ax.finallimits, ax.limits),
        shiftax_limits = obs_getset(shiftax.finallimits, shiftax.limits),
    )
    get_cur_state() = modify(o -> o[], stateobs, RecursiveOfType(Observable))
    set_cur_state(vals) = modify((o, v) -> o[] = v, stateobs, RecursiveOfType(Observable), vals)
    statedefault = get_cur_state()
    on(Button(fig[0,:][1,3], label="Reset").clicks) do _
        set_cur_state(statedefault)
    end
    state_db_file = joinpath(@get_scratch!("cache"), "cache.db")
    state_db = SQLite.DB(state_db_file)
    state_dict = SQLDictionary{@NamedTuple{f1::String, f2::String}, @NamedTuple{ui::Any}}(state_db, "ui_state")
    stored_state = get(state_dict, statekey, nothing)
    if !isnothing(stored_state)
        if typeof(stored_state.ui) == typeof(statedefault)
            set_cur_state(stored_state.ui)
        else
            @warn "Previous UI state is stored, but types do not match. Resetting to default."
        end
    end
    stateobs_trigger = lift(Returns(nothing), stateobs...)
    on(Makie.Observables.throttle(2, stateobs_trigger)) do _
        set!(state_dict, statekey, (;ui=get_cur_state()))
    end

    if display
        MakieExtra.show_gl_icon_in_dock()
        Base.display(fig)
    else
        return fig
    end
end

function obs_getset(obs_get, obs_set)
    obs = typeof(obs_get)(obs_get[])
    selfupdating = Base.Lockable(Ref(false))
    on(obs_get) do v
        lock(selfupdating) do selfupdating
            selfupdating[] = true
            obs[] = v
            selfupdating[] = false
        end
    end
    on(obs) do v
        lock(selfupdating) do selfupdating
            if !selfupdating[]
                obs_set[] = v
            end
        end
    end
    return obs
end


end

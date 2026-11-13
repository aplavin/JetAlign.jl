using TestItems
using TestItemRunner
@run_package_tests


@testitem "align_strings" begin
    using JetAlign: align_strings, difference_runs

    # Test exact match
    @test align_strings("hello", "hello") == [(:match, 'h', 'h'), (:match, 'e', 'e'), (:match, 'l', 'l'), (:match, 'l', 'l'), (:match, 'o', 'o')]

    # Test substitution
    @test align_strings("hello", "hallo") == [(:match, 'h', 'h'), (:substitute, 'e', 'a'), (:match, 'l', 'l'), (:match, 'l', 'l'), (:match, 'o', 'o')]

    # Test insertion
    @test align_strings("cat", "cart") == [(:match, 'c', 'c'), (:match, 'a', 'a'), (:insert, nothing, 'r'), (:match, 't', 't')]

    # Test deletion
    @test align_strings("testing", "test") == [(:match, 't', 't'), (:match, 'e', 'e'), (:match, 's', 's'), (:match, 't', 't'), (:delete, 'i', nothing), (:delete, 'n', nothing), (:delete, 'g', nothing)]

    # Test empty strings
    @test align_strings("", "") == []

    # Test one empty string
    @test align_strings("abc", "") == [(:delete, 'a', nothing), (:delete, 'b', nothing), (:delete, 'c', nothing)]

    @test align_strings("", "xyz") == [(:insert, nothing, 'x'), (:insert, nothing, 'y'), (:insert, nothing, 'z')]

    # Test complex case with multiple operations
    @test align_strings("kitten", "sitting") == [(:substitute, 'k', 's'), (:match, 'i', 'i'), (:match, 't', 't'), (:match, 't', 't'), (:substitute, 'e', 'i'), (:match, 'n', 'n'), (:insert, nothing, 'g')]

    # Test Cyrillic characters
    @test align_strings("привет", "привет") == [(:match, 'п', 'п'), (:match, 'р', 'р'), (:match, 'и', 'и'), (:match, 'в', 'в'), (:match, 'е', 'е'), (:match, 'т', 'т')]

    @test align_strings("кот", "кит") == [(:match, 'к', 'к'), (:substitute, 'о', 'и'), (:match, 'т', 'т')]

    @test align_strings("мир", "мира") == [(:match, 'м', 'м'), (:match, 'и', 'и'), (:match, 'р', 'р'), (:insert, nothing, 'а')]

    # Test mixed scripts
    @test align_strings("hello мир", "hello мира") == [(:match, 'h', 'h'), (:match, 'e', 'e'), (:match, 'l', 'l'), (:match, 'l', 'l'), (:match, 'o', 'o'), (:match, ' ', ' '), (:match, 'м', 'м'), (:match, 'и', 'и'), (:match, 'р', 'р'), (:insert, nothing, 'а')]

    # Test other Unicode characters (emoji, mathematical symbols)
    @test align_strings("😀😃", "😀😄") == [(:match, '😀', '😀'), (:substitute, '😃', '😄')]

    @test align_strings("α + β", "α - β") == [(:match, 'α', 'α'), (:match, ' ', ' '), (:substitute, '+', '-'), (:match, ' ', ' '), (:match, 'β', 'β')]

    # Test Chinese/Japanese characters
    @test align_strings("你好", "你好") == [(:match, '你', '你'), (:match, '好', '好')]

    @test align_strings("こんにちは", "こんばんは") == [(:match, 'こ', 'こ'), (:match, 'ん', 'ん'), (:substitute, 'に', 'ば'), (:substitute, 'ち', 'ん'), (:match, 'は', 'は')]

    # runs of differing characters, as displayed for the two file names
    @test difference_runs("J0433_S_map.fits", "J0433_X_map.fits") == (
        [(text="J0433_", differs=false), (text="S", differs=true), (text="_map.fits", differs=false)],
        [(text="J0433_", differs=false), (text="X", differs=true), (text="_map.fits", differs=false)])
    @test difference_runs("cat", "cart") == ([(text="cat", differs=false)], [(text="ca", differs=false), (text="r", differs=true), (text="t", differs=false)])
    @test difference_runs("", "ab") == ([], [(text="ab", differs=true)])
end

@testitem "spectral index" begin
    using JetAlign: spectral_index
    using JetAlign.Unitful

    a = [1.0 2.0; 4.0 0.03]
    b = a .* 2.0^-0.7                       # S ∝ ν^-0.7 between 1 and 2 GHz
    si = spectral_index(a, b, (1u"GHz", 2u"GHz"); dynrange=100)
    @test si[1:3] ≈ fill(-0.7, 3) rtol=1e-12
    @test isnan(si[4])                      # below 4/100 in `a`
    # either image too faint masks the pixel
    @test isnan.(spectral_index([1.0, 1.0], [1.0, 1e-3], (1u"GHz", 2u"GHz"); dynrange=100)) == [false, true]
end

@testitem "cc_to_image" begin
    using JetAlign: cc_to_image, skygrid
    using JetAlign.InterferometricModels, JetAlign.Unitful, JetAlign.UnitfulAngles, JetAlign.StaticArrays, JetAlign.AxisKeysExtra

    model = MultiComponentModel([
        InterferometricModels.Point(1.0, SVector(0.0, 0.0)u"mas"),
        InterferometricModels.Point(0.5, SVector(2.0, -1.0)u"mas"),
        InterferometricModels.Point(7.0, SVector(9.0, 0.0)u"mas"),   # outside the grid: dropped
    ])
    radecs = skygrid(5, 101)                 # 0.1 mas pixels
    img = cc_to_image(model; σ=0.3u"mas", radecs)
    # peaks equal the fluxes (Jy/beam), at the component positions; the other component adds 0.5·exp(-5/0.18)
    @test img(ra=0u"mas", dec=0u"mas") ≈ 1.0 rtol=1e-12
    @test img(ra=2u"mas", dec=-1u"mas") ≈ 0.5 rtol=1e-11
    @test argmax(img) == CartesianIndex(51, 51)
    # total: Σ 2πσ²F / pixel area
    @test sum(img) ≈ 2π * 0.09 * 1.5 / 0.01 rtol=1e-9
    @test maximum(img(ra=5u"mas")) < 1e-20           # nothing piled up at the edge

    shifted = cc_to_image(model; σ=0.3u"mas", radecs, shift=SVector(0.5, -0.2)u"mas")
    @test argmax(shifted) == CartesianIndex(56, 49)
    @test maximum(shifted) ≈ 1.0 rtol=1e-12
end

@testitem "x-corr peak" begin
    using JetAlign: xcorr_ccs, skygrid
    using JetAlign.InterferometricModels, JetAlign.Unitful, JetAlign.UnitfulAngles, JetAlign.StaticArrays, JetAlign.AxisKeysExtra

    pts = [(1.0, SVector(0.0, 0.0)), (0.5, SVector(1.0, 0.0)), (0.3, SVector(0.0, 1.2))]
    a = MultiComponentModel([InterferometricModels.Point(f, p * u"mas") for (f, p) in pts])
    b = MultiComponentModel([InterferometricModels.Point(f, (p + SVector(0.4, -0.3)) * u"mas") for (f, p) in pts])
    xc = xcorr_ccs(b, a; radecs=skygrid(2, 41))
    i, j = Tuple(argmax(xc))
    @test axiskeys(xc, :ra)[i] ≈ 0.4u"mas"
    @test axiskeys(xc, :dec)[j] ≈ -0.3u"mas"
    @test maximum(xc) ≈ 1.0^2 + 0.5^2 + 0.3^2 rtol=1e-6
end

@testitem "beam polygon" begin
    using JetAlign: beampoly
    using JetAlign.InterferometricModels, JetAlign.Unitful, JetAlign.UnitfulAngles

    b = beam(EllipticGaussian; σ_major=1.0u"mas", ratio_minor_major=0.5, pa_major=π/6)
    pts = beampoly(b; n=5)                  # θ = 0, π/2, π, 3π/2, 2π
    a = 2.3548200450309493 / 2              # FWHM/2 of σ = 1
    @test pts[1] ≈ a * [sin(π/6), cos(π/6)] * u"mas"         # along the position angle, from north to east
    @test pts[2] ≈ a / 2 * [cos(π/6), -sin(π/6)] * u"mas"    # minor axis
    @test pts[3] ≈ -pts[1]
    @test pts[5] ≈ pts[1]
    @test map(p -> hypot(p...), beampoly(beam(CircularGaussian; σ=1.0u"mas"))) ≈ fill(a, 65)u"mas"
end

@testitem "state persistence" begin
    using JetAlign: UIState, open_store, restore, save!
    using JetAlign.AccessorsExtra, JetAlign.IntervalSets, JetAlign.StaticArrays

    lims = (ra=-1.0..1.0, dec=-2.0..2.0)
    default = UIState(; fwhm=1.0, maskradius=1.5, halffov=10.0, shiftlimits=lims, skylimits=lims)
    changed = setproperties(default; npix=256, shift=SVector(0.3, -0.1), sicolorrange=(-1.0, 1.0), skylimits=(ra=-3.0..2.5, dec=-1.0..1.5))
    mktempdir() do dir
        db = joinpath(dir, "state.db")
        files = ["a.fits", "b.fits"]
        @test restore(open_store(db, files), default) === default
        save!(open_store(db, files), changed)
        @test restore(open_store(db, files), default) === changed
        @test restore(open_store(db, reverse(files)), default) === default
        # a state stored by another version of UIState is not restored
        store = open_store(db, files)
        JetAlign.SQLCollections.set!(store.dict, store.key, (; ui=(; fwhm=1.0)))
        @test (@test_logs (:warn,) restore(store, default)) === default
    end
    @test restore(nothing, default) === default
end

@testitem "shifts file" begin
    using JetAlign: save_shifts, load_shifts, App
    using JetAlign.StaticArrays, JetAlign.VLBIFiles

    shifts = [SVector(0.5, -0.25), SVector(-1.0, 2.0)]
    mktempdir() do dir
        file = joinpath(dir, "shifts.toml")
        save_shifts(file, ["a.fits", "b.fits"], shifts)
        @test load_shifts(file, ["a.fits", "b.fits"]) == shifts
        @test_throws "not for" load_shifts(file, ["b.fits", "a.fits"])
        @test_throws "not for" load_shifts(file, ["a.fits", "c.fits"])

        map_fits = joinpath(pkgdir(VLBIFiles), "test", "data", "map.fits")
        @test_throws "not for" App([map_fits, map_fits]; shiftsfile=file, persiststate=false)
        save_shifts(file, [map_fits, map_fits], shifts)
        @test App([map_fits, map_fits]; shiftsfile=file, persiststate=false).shifts == shifts
    end
end

# Drives the real window under ImGuiTestEngine on the VLBIFiles map.fits pair: every slider, the shift
# DragPoint, A/D+click on multishift points, Save shifts and Reset. Needs a display.
@testitem "windowed app" begin
    using JetAlign: App, run!, maps_inputs, xcorr_inputs, shifts_path, SingleSI, MultiSI, NPIX_CHOICES
    using JetAlign.StaticArrays, JetAlign.VLBIFiles, JetAlign.IntervalSets
    import ImGuiTestEngine as te
    using ImGuiTestEngine: @register_test, @imcheck
    import CImGui as ig
    import ImPlot
    import TOML

    map_fits = joinpath(pkgdir(VLBIFiles), "test", "data", "map.fits")
    app = App([map_fits, map_fits]; persiststate=false)

    engine = te.CreateContext(; exit_on_completion=true, show_test_window=false)
    engine.IO.ConfigWatchdogWarning = 3600f0       # first frames compile for a while
    engine.IO.ConfigWatchdogKillTest = 3600f0

    # yield frames until `pred()`, or the frame budget runs out; returns `pred()`
    function yield_until(pred; budget=20000)
        for _ in 1:budget
            pred() && return true
            te.Yield()
        end
        pred()
    end
    settled(job, inputs) = !isnothing(job.published) && job.published.inputs == inputs && isnothing(job.task)
    jobs_settled() = settled(app.jobs.maps, maps_inputs(app)) && settled(app.jobs.xcorr, xcorr_inputs(app))

    # screen position of data point `(x, y)` in the shift plot: the leftmost of the app's plots
    function shiftplot_pixel(x, y)
        pool = unsafe_load(ImPlot.GetCurrentContext()).Plots.Buf
        plot = argmin(p -> p.FrameRect.Min.x, [unsafe_load(pool.Data, i) for i in 1:pool.Size])
        ax, ay = plot.Axes[Int(ImPlot.ImAxis_X1) + 1], plot.Axes[Int(ImPlot.ImAxis_Y1) + 1]
        ig.ImVec2(ax.PixelMin + ax.ScaleToPixel * (x - ax.ScaleMin), ay.PixelMin + ay.ScaleToPixel * (y - ay.ScaleMin))
    end
    function hotkey_click(key, p)
        te.MouseMoveToPos(shiftplot_pixel(p...))
        te.KeyDown(key); te.Yield(2)
        te.MouseClick(ig.ImGuiMouseButton_Left)
        te.KeyUp(key); te.Yield(2)
    end

    results = Dict{Symbol,Any}()
    @register_test(engine, "JetAlign", "e2e") do
        te.Yield(10)
        @imcheck yield_until(jobs_settled)
        te.SetRef("JetAlign")

        # every slider
        te.ItemInputValue("Beam FWHM", 3.0)
        te.ItemInputValue("Dynamic range", 1000.0)
        te.ItemInputValue("npix", 6)
        te.ItemInputValue("FoV", 20.0)
        te.ItemInputValue("X-corr mask radius", 1.0)
        te.Yield(2)
        results[:sliders] = (app.ui.fwhm, app.ui.dynrange, app.ui.npix, app.ui.halffov, app.ui.maskradius)
        results[:recomputed] = yield_until(jobs_settled)
        results[:maps] = app.jobs.maps.published.result

        # drag the shift point from the origin
        te.MouseMoveToPos(shiftplot_pixel(0, 0))
        te.MouseDown(ig.ImGuiMouseButton_Left)
        te.MouseMoveToPos(shiftplot_pixel(1.0, 0.5))
        te.MouseUp(ig.ImGuiMouseButton_Left)
        te.Yield(2)
        results[:shift] = (app.ui.shift, app.simode)

        # A+click adds multishift points, D+click deletes the nearest one
        hotkey_click(ig.ImGuiKey_A, (0.5, -0.5))
        hotkey_click(ig.ImGuiKey_A, (-0.5, 0.5))
        results[:added] = (copy(app.shifts), app.simode, app.ui.shift)
        hotkey_click(ig.ImGuiKey_D, (0.6, -0.4))
        results[:deleted] = (copy(app.shifts), app.simode)
        @imcheck yield_until(jobs_settled)
        results[:si] = app.jobs.maps.published.result.si
        hotkey_click(ig.ImGuiKey_D, (0, 0))
        results[:alldeleted] = (copy(app.shifts), app.simode)
        hotkey_click(ig.ImGuiKey_A, (-0.5, 0.5))

        te.ItemClick("Save shifts")
        results[:saved] = TOML.parsefile(shifts_path(app))

        te.ItemClick("Reset")
        te.Yield(5)
        results[:reset] = (app.ui, app.simode)
    end

    run!(app; engine)
    rm(shifts_path(app); force=true)

    summary = te.GetResultSummary(engine)
    @test summary.CountTested == 1
    @test summary.CountSuccess == 1

    @test collect(results[:sliders]) ≈ [3.0, 1000.0, NPIX_CHOICES[7], 20.0, 1.0] rtol=1e-6
    @test results[:recomputed]
    @test size(results[:maps].imgs[1]) == (256, 256)
    @test collect(endpoints(results[:maps].extent[1])) == [-20, 20]
    @test results[:maps].peaks[1] ≈ results[:maps].peaks[2]

    let (shift, mode) = results[:shift]
        @test shift ≈ [1.0, 0.5] atol=0.05
        @test mode === SingleSI()
    end
    let (shifts, mode, shift) = results[:added]
        @test length(shifts) == 2
        @test shifts[1] ≈ [0.5, -0.5] atol=0.05
        @test shifts[2] ≈ [-0.5, 0.5] atol=0.05
        @test mode === MultiSI()
        @test shift ≈ results[:shift][1]               # the clicks did not move the shift point
    end
    let (shifts, mode) = results[:deleted]
        @test shifts == results[:added][1][2:2]
        @test mode === MultiSI()
    end
    @test eltype(results[:si]) <: ClosedInterval
    @test results[:alldeleted] == (SVector{2,Float64}[], SingleSI())

    @test results[:saved]["files"] == [map_fits, map_fits]
    @test only(results[:saved]["shifts"]) ≈ [-0.5, 0.5] atol=0.05

    let (ui, mode) = results[:reset]
        d = app.defaults
        @test (ui.fwhm, ui.dynrange, ui.maskradius, ui.npix, ui.halffov, ui.shift, ui.sicolorrange) ==
              (d.fwhm, d.dynrange, d.maskradius, d.npix, d.halffov, d.shift, d.sicolorrange)
        @test collect(endpoints(ui.shiftlimits.dec)) ≈ collect(endpoints(d.shiftlimits.dec))   # equal aspect: ImPlot fits RA to Dec
        @test collect(endpoints(ui.skylimits.dec)) ≈ collect(endpoints(d.skylimits.dec))
        @test mode === SingleSI()
    end

    te.DestroyContext(engine)
end

@testitem "_" begin
    import Aqua
    Aqua.test_all(JetAlign)

    import CompatHelperLocal as CHL
    CHL.@check()
end

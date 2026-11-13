imvec4(c::Colorant) = (c = RGBA{Float32}(c); ig.ImVec4(red(c), green(c), blue(c), alpha(c)))
imu32(c::Colorant) = ig.GetColorU32(imvec4(c))

const PLOT_FLAGS = ImPlot.ImPlotFlags_Equal | ImPlot.ImPlotFlags_NoLegend
const NOFIT = Int(ImPlot.ImPlotItemFlags_NoFit)
const SHIFT_COLOR = colorant"white"
const SHIFTS_COLOR = colorant"#0072b2"

function init_context!()
    ig.set_backend(:GlfwOpenGL3)
    ctx = ig.CreateContext()
    pctx = ImPlot.CreateContext()
    ImPlot.SetImGuiContext(ctx)
    ig.GetIO().IniFilename = C_NULL      # nothing to persist: one fixed window
    ig.AddFontFromFileTTF(ImGuiThemes.FONTS[(:JuliaMono, :Regular)])
    ImGuiThemes.apply_theme!("Solarized Light")
    ctx, pctx
end

"Open the app window; returns when it is closed. `engine`: an ImGuiTestEngine context driving it."
function run!(app::App; engine=nothing)
    ctx, pctx = init_context!()
    ig.render(ctx; engine, spawn=:interactive, window_title="JetAlign.jl", window_size=(1400, 650),
              on_exit=() -> (persist!(app; interval=0); ImPlot.DestroyContext(pctx))) do
        frame!(app)
    end
    nothing
end

function frame!(app::App)
    maps = tick!(inputs -> compute_maps(app.images, inputs), app.jobs.maps, maps_inputs(app))
    xcorr = tick!(inputs -> compute_xcorr(app.images, app.xcorr_halfwidth, inputs), app.jobs.xcorr, xcorr_inputs(app))
    vp = ig.GetMainViewport()
    ig.SetNextWindowPos(unsafe_load(vp.WorkPos))
    ig.SetNextWindowSize(unsafe_load(vp.WorkSize))
    flags = ig.ImGuiWindowFlags_NoDecoration | ig.ImGuiWindowFlags_NoMove | ig.ImGuiWindowFlags_NoSavedSettings |
            ig.ImGuiWindowFlags_NoBringToFrontOnFocus
    if ig.Begin("JetAlign", C_NULL, flags)
        controls!(app)
        plots!(app, maps, xcorr)
    end
    ig.End()
    persist!(app; interval=2)
end


# ---------------- controls ----------------

function controls!(app::App)
    ig.BeginGroup()
    for (color, runs) in zip(app.colors, difference_runs(map(im -> im.file, app.images)...))
        for (i, r) in enumerate(runs)
            i > 1 && ig.SameLine(0, 0)
            r.differs && ig.PushStyleColor(ig.ImGuiCol_Text, imvec4(color))
            ig.TextUnformatted(r.text)
            r.differs && ig.PopStyleColor()
        end
    end
    ig.Button("Reset") && reset!(app)
    ig.SameLine()
    ig.Button("Save shifts") && save_shifts!(app)
    ig.EndGroup()

    ig.SameLine(0, 40)
    ig.BeginGroup()
    ig.PushItemWidth(200)
    ui, σ0 = app.ui, fwhm_to_σ(app.defaults.fwhm)
    ui = @set ui.fwhm = slider("Beam FWHM", ui.fwhm, 0.05, 10σ0, "%.1f mas")
    ig.SameLine()
    ui = @set ui.dynrange = slider("Dynamic range", ui.dynrange, 10, 1e6, "%.2g", ig.ImGuiSliderFlags_Logarithmic)
    ui = @set ui.npix = choice_slider("npix", NPIX_CHOICES, ui.npix)
    ig.SameLine()
    ui = @set ui.halffov = slider("FoV", ui.halffov, 0.5, 200, "%.1f mas", ig.ImGuiSliderFlags_Logarithmic)
    ui = @set ui.maskradius = slider("X-corr mask radius", ui.maskradius, 0, 25σ0, "%.1f mas")
    ig.PopItemWidth()
    ig.EndGroup()
    app.ui = ui
end

"Float slider showing `x`; returns the new value, or `x` itself when untouched."
function slider(label, x, lo, hi, fmt, flags=0)
    v = Ref(Cfloat(x))
    ig.SliderFloat(label, v, lo, hi, fmt, flags) ? Float64(v[]) : x
end

"Slider over the index of `x` in `values`; returns the selected value."
function choice_slider(label, values, x)
    i = Ref(Cint(findfirst(==(x), values) - 1))
    ig.SliderInt(label, i, 0, length(values) - 1, string(x))
    values[i[] + 1]
end


# ---------------- plots ----------------

function plots!(app::App, maps, xcorr)
    ig.BeginTable("plots", 4) || return
    foreach(("shift", "I", "SI")) do name
        ig.TableSetupColumn(name, ig.ImGuiTableColumnFlags_WidthStretch)
    end
    ig.TableSetupColumn("colorbar", ig.ImGuiTableColumnFlags_WidthFixed)
    ig.TableNextRow()
    h = ig.GetContentRegionAvail().y
    ig.TableNextColumn()
    shift_plot!(app, xcorr, h)

    ui = app.ui
    sky = limitrefs(ui.skylimits)
    cr = Ref(ui.sicolorrange)
    GC.@preserve sky cr begin
        ig.TableNextColumn()
        hovered_stokes = stokes_plot!(app, maps, sky, h)
        ig.TableNextColumn()
        hovered_si = si_plot!(app, maps, sky, h)
        ig.TableNextColumn()
        colorbar!("##sibar", :turbo, cr, identity, h)
        app.cursor = something(hovered_stokes, hovered_si, SVector(NaN, NaN))
        app.ui = setproperties(ui; skylimits=readlimits(sky), sicolorrange=cr[])
    end
    ig.EndTable()
end

limitrefs(l::Limits) = map(iv -> (Ref(leftendpoint(iv)), Ref(rightendpoint(iv))), l)
readlimits(refs) = map(((lo, hi),) -> lo[]..hi[], refs)
link_axes!(refs) = (ImPlot.SetupAxisLinks(ImPlot.ImAxis_X1, refs.ra...); ImPlot.SetupAxisLinks(ImPlot.ImAxis_Y1, refs.dec...))

"Mouse position in the current plot when it is hovered, else nothing."
plot_mouse() = ImPlot.IsPlotHovered() ? SVector(ImPlot.GetPlotMousePos().x, ImPlot.GetPlotMousePos().y) : nothing

"Data position at fractions `(fx, fy)` of the current plot area, from its bottom-left corner on screen."
function plot_point(fx, fy)
    pos, size = ImPlot.GetPlotPos(), ImPlot.GetPlotSize()
    p = ImPlot.PixelsToPlot(ig.ImVec2(pos.x + fx * size.x, pos.y + (1 - fy) * size.y))
    SVector(p.x, p.y)
end

"Lines of `(text, color::ImU32)` in the top-left or top-right corner of the current plot."
function plot_note!(lines; right=false)
    pos, size = ImPlot.GetPlotPos(), ImPlot.GetPlotSize()
    pad = 5
    foreach(enumerate(lines)) do (i, (text, color))
        x = right ? pos.x + size.x - pad - ig.CalcTextSize(text).x : pos.x + pad
        y = pos.y + pad + (i - 1) * ig.GetTextLineHeight()
        ig.AddText(ImPlot.GetPlotDrawList(), ig.ImVec2(x, y), color, text)
    end
end

"Closed outline through `points` (angular `SVector`s), offset by `center` (in mas)."
function outline!(label, points, center, color; weight=1)
    xy = map(p -> ustrip.(u"mas", p) + center, points)
    ImPlot.PlotLine(label, first.(xy), last.(xy); spec=ImPlot.ImPlotSpec(; LineColor=imvec4(color), LineWeight=weight, Flags=NOFIT))
end

curbeam(app::App) = beam(CircularGaussian, σ=fwhm_to_σ(app.ui.fwhm) * u"mas")

function shift_plot!(app::App, xcorr, h)
    # A/D + click edit the multishift points: they own the left button so it neither pans nor drags
    hotkey = ig.IsKeyDown(ig.ImGuiKey_A) || ig.IsKeyDown(ig.ImGuiKey_D)
    owner = ig.GetID("shift hotkeys")
    hotkey && ig.SetKeyOwner(ig.ImGuiKey_MouseLeft, owner)
    lims = limitrefs(app.ui.shiftlimits)
    ImPlot.PushStyleColor(ImPlot.ImPlotCol_PlotBg, imvec4(colorant"black"))
    GC.@preserve lims if ImPlot.BeginPlot("##shift", "RA (mas)", "Dec (mas)", ig.ImVec2(-1, h);
                                         flags=PLOT_FLAGS, x_flags=ImPlot.ImPlotAxisFlags_Invert)
        link_axes!(lims)
        if !isnothing(xcorr)
            (; extent, img, colorrange) = xcorr.result
            image!("##xcorr", extent..., img; colormap=:turbo, colorrange, interpolate=true, flags=NOFIT)
        end
        white = imvec4(colorant"white")
        ImPlot.PlotInfLines("##zero", [0.0]; spec=ImPlot.ImPlotSpec(; LineColor=white, Flags=NOFIT))
        ImPlot.PlotInfLines("##zero", [0.0]; spec=ImPlot.ImPlotSpec(; LineColor=white, Flags=NOFIT | Int(ImPlot.ImPlotInfLinesFlags_Horizontal)))
        outline!("##beam", beampoly(curbeam(app)), plot_point(0.1, 0.1), RGBA(1, 1, 1, 0.3))
        plot_note!([(text, imu32(colorant"white")) for text in
                    ("Color: x-corr", "Drag: image shift", "A+click: add point", "D+click: delete point")]; right=true)

        dragflags = hotkey ? ImPlot.ImPlotDragToolFlags_NoInputs : ImPlot.ImPlotDragToolFlags_None
        drag_point(0, app.ui.shift, SHIFT_COLOR, dragflags) do p
            set_shift!(app, p)
        end
        foreach(enumerate(app.shifts)) do (i, s)
            drag_point(i, s, SHIFTS_COLOR, dragflags) do p
                shifts = app.shifts
                set_shifts!(app, @set shifts[i] = p)
            end
        end
        if hotkey && ImPlot.IsPlotHovered() && ig.IsMouseClicked(ig.ImGuiMouseButton_Left, ig.ImGuiInputFlags_None, owner)
            p = plot_mouse()
            if ig.IsKeyDown(ig.ImGuiKey_A)
                set_shifts!(app, [app.shifts; [p]])
            elseif !isempty(app.shifts)
                set_shifts!(app, deleteat!(copy(app.shifts), argmin(i -> norm(app.shifts[i] - p), eachindex(app.shifts))))
            end
        end
        ImPlot.EndPlot()
        app.ui = setproperties(app.ui; shiftlimits=readlimits(lims))
    end
    ImPlot.PopStyleColor()
end

"`ImPlot.DragPoint` at `p`; calls `onmove(new position)` when dragged."
function drag_point(onmove, id, p, color, flags)
    x, y = Ref(p[1]), Ref(p[2])
    ImPlot.DragPoint(id, x, y, imvec4(color), 5, flags) && onmove(SVector(x[], y[]))
end

function stokes_plot!(app::App, maps, sky, h)
    ImPlot.BeginPlot("##stokes", "", "", ig.ImVec2(-1, h); flags=PLOT_FLAGS, x_flags=ImPlot.ImPlotAxisFlags_Invert) || return nothing
    link_axes!(sky)
    if !isnothing(maps)
        (; extent, imgs, levels) = maps.result
        foreach(1:2) do i
            contour!("##I$i", extent..., imgs[i]; levels=levels[i], color=app.colors[i], flags=NOFIT)
        end
    end
    r = app.ui.maskradius
    mask = ellipse(r, r, 0)
    ImPlot.PlotPolygon("##mask", first.(mask), last.(mask);
        spec=ImPlot.ImPlotSpec(; FillColor=imvec4(colorant"black"), FillAlpha=0.12, LineColor=imvec4(colorant"black"), LineWeight=0.5, Flags=NOFIT))
    corner = plot_point(0.1, 0.1)
    outline!("##beam", beampoly(curbeam(app)), corner, colorant"black")
    foreach(enumerate(app.images)) do (i, im)
        outline!("##beam$i", beampoly(im.beam), corner, RGBA(app.colors[i], 0.3))
    end
    plot_cursor!(app)
    plot_note!([("Stokes I", ig.GetColorU32(ig.ImGuiCol_Text))])
    if !isnothing(maps)
        plot_note!([(f"{im.freq:.1f}: {peak:.2f} Jy/bm", imu32(app.colors[i]))
                    for (i, (im, peak)) in enumerate(zip(app.images, maps.result.peaks))]; right=true)
    end
    hovered = plot_mouse()
    ImPlot.EndPlot()
    hovered
end

function si_plot!(app::App, maps, sky, h)
    ImPlot.BeginPlot("##si", "", "", ig.ImVec2(-1, h); flags=PLOT_FLAGS, x_flags=ImPlot.ImPlotAxisFlags_Invert) || return nothing
    link_axes!(sky)
    isnothing(maps) || si_image!(app, maps.result.extent, maps.result.si)
    plot_cursor!(app)
    plot_note!([("Spectral index", ig.GetColorU32(ig.ImGuiCol_Text))]; right=true)
    hovered = plot_mouse()
    ImPlot.EndPlot()
    hovered
end

si_image!(app::App, extent, si::AbstractMatrix{<:Real}) =
    image!("##si", extent..., si; colormap=:turbo, colorrange=app.ui.sicolorrange, interpolate=true, flags=NOFIT)
function si_image!(app::App, extent, si::AbstractMatrix{<:ClosedInterval})
    cr = app.ui.sicolorrange
    colors = cached(() -> map(int -> si_color(cr, int), si), (si, cr), app.sicolors)
    image!("##si", extent..., colors; interpolate=true, flags=NOFIT)
end

plot_cursor!(app::App) = all(isfinite, app.cursor) &&
    ImPlot.PlotScatter("##cursor", [app.cursor[1]], [app.cursor[2]]; spec=ImPlot.ImPlotSpec(;
        Marker=ImPlot.ImPlotMarker_Cross, MarkerSize=8, LineWeight=2, MarkerLineColor=imvec4(colorant"black"), Flags=NOFIT))

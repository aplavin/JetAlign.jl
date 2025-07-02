module JetAlign


using DataManipulation
# using StructArrays
using StaticArrays
using AccessorsExtra
using InterferometricModels
using VLBIFiles; using VLBIFiles: coords, flux
using AxisKeysExtra
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


round_(x::Number, precision::Number) = round(Integer, x / precision) * precision


mask_ccs(mod::MultiComponentModel, radius) = @modify(components(mod)) do ccs
    filter(ccs) do cc
        norm(coords(cc)) > radius
    end
end

module CorrelationKind
struct Raw end
struct Normalized end
end

function aggregate_ccs_by_pixel(ccs; step)
    @p @modify(components(ccs)) do __
        group_vg(round_.(coords(_), step))
        map(InterferometricModels.Point(sum(flux, _), key(_)))
    end
end

function spectral_index(a::KeyedArray, b::KeyedArray, freqs; threshold=0, threshold_dynrange=Inf)
    @assert length(freqs) == 2
    shifted_imgs = a, b
    log_freq_ratio = log(freqs[2] / freqs[1])
    si = similar(shifted_imgs[2])
    si .= log.(clamp.(shifted_imgs[2] ./ shifted_imgs[1], 0, Inf)) ./ log_freq_ratio
    # shifted_imgs[2]
    for sh_img in shifted_imgs
        peak = maximum(sh_img)
        si[sh_img .< max(threshold, peak / threshold_dynrange)] .= NaN
    end
    return si
end

function smooth!(A::Array; σ_pix::Real)
    pad = ceil(Int, 1*σ_pix)
    A[begin:begin+pad, :] .= 0
    A[end-pad:end, :] .= 0
    A[:, begin:begin+pad] .= 0
    A[:, end-pad:end] .= 0

    kernel = Kernel.gaussian((σ_pix,σ_pix), (8*ceil(Int, σ_pix) + 1, 8*ceil(Int, σ_pix) + 1))

    tmp = map(Returns(0.), A)
    tmp[end ÷ 2, end ÷ 2] = 1
    imfilter!(tmp, tmp, kernel)	

    imfilter!(A, A, kernel)
    A ./= maximum(tmp)
    # imfilter!(A, A, IF.KernelFactors.IIRGaussian((σ_pix, σ_pix)))

    return A
end

function smooth!(A::KeyedArray; σ::Number)
    pixstep = @p axiskeys(A) map(step(_) |> abs) uniqueonly
    σ_pix = σ / pixstep
    smooth!(A.data.data; σ_pix)
    return A
end

smooth(A::AbstractArray; kwargs...) = smooth!(copy(A); kwargs...)

function unsmoothed_xcorr_ccs(::CorrelationKind.Raw, ccs_a, ccs_b; radecs)
    pixstep = @p axiskeys(radecs) map(step(_) |> abs) uniqueonly
    xcorr = map(Returns(0.), radecs)
    ccs_a_ = aggregate_ccs_by_pixel(ccs_a; step=pixstep/5)
    ccs_b_ = aggregate_ccs_by_pixel(ccs_b; step=pixstep/5)
    for cc_a in components(ccs_a_)
        for cc_b in components(ccs_b_)
            diff = coords(cc_a) - coords(cc_b)
            xcorr(Near(diff[1]), Near(diff[2]), :) .+= flux(cc_a) * flux(cc_b)
        end
    end
    return xcorr
end

function xcorr_ccs(kind::CorrelationKind.Raw, ccs_a, ccs_b; radecs, σ)
    smooth!(unsmoothed_xcorr_ccs(kind, ccs_a, ccs_b; radecs); σ)
end

function cc_to_image(cc; σ, radecs, shift=SVector(0u"mas", 0u"mas"))
    img = map(Returns(0.), radecs)
    for c in components(cc)
        coo = coords(c) + shift
        img(Near(coo[1]), Near(coo[2]), :) .+= c.flux
    end
    return smooth!(img; σ)
end

function spectral_index(a::MultiComponentModel, b::MultiComponentModel, freqs; shift, radecs, σ, kwargs...)
    ccs_ab = (a, b)
    shifted_imgs = (cc_to_image(ccs_ab[1]; σ, radecs, shift), cc_to_image(ccs_ab[2]; σ, radecs))
    spectral_index(shifted_imgs..., freqs; kwargs...)
end

function highlight_string_differences(strs; colors)
	if allequal(length, strs)
		equalchars = @p let
			strs
			map() do str
				collect(str) .== collect(first(strs))
			end
			stack
			eachslice(dims=1)
			map(all)
		end
		map(strs, colors) do str, color
			parts = map(collect(str), equalchars) do ch, isequal
				kws = isequal ? (;) : (;color)
				rich(string(ch); kws...)
			end
			rich(parts...)
		end
	else
		neqchars_start = @p let
			0:minimum(length, strs)
			filter() do i
				allequal(s -> (@view s[begin:i]), strs)
			end
			maximum
		end
		neqchars_end = @p let
			0:minimum(length, strs)
			filter() do i
				allequal(s -> (@view s[end-i+1:end]), strs)
			end
			minimum
		end
		map(strs, colors) do str, color
			rich(
				str[1:neqchars_start],
				rich(str[neqchars_start+1:end-neqchars_end]; color),
				str[end-neqchars_end+1:end],
			)
		end
	end
end

function imalign_interactive(files)
    @assert length(files) == 2
    files = map(expanduser, files)
    @assert all(isfile, files)

    images_with_data = @p map(files) do f
        fimg = VLBI.load(f, read_data=false)
        mod = VLBI.load(MultiComponentModel, f)
        mod = @modify(ustrip, mod |> RecursiveOfType(ModelComponent) |> flux)
        (;freq=VLBI.frequency(fimg), mod, beam=VLBI.load(Beam, f))
    end sort(by=_.freq)

    colors = [:red, :blue]
    radecs = @p range(-8u"mas", 8u"mas", length=512) grid(ra=reverse(__), dec=__)

    with_theme(Axis=(;width=350, height=350)) do
        fig = Figure()

        σobs = Observable(0.5)
        curbeam = @lift beam(CircularGaussian, σ=$σobs)

        shiftobs = Observable(SVector(0., 0.))
        img_a = @lift cc_to_image(images_with_data[1].mod; σ=$σobs*u"mas", radecs, shift=SVector($shiftobs)u"mas")
        img_b = @lift cc_to_image(images_with_data[2].mod; σ=$σobs*u"mas", radecs)
        SI = @lift spectral_index($img_a, $img_b, map(:freq, images_with_data); threshold=0.3e-3)

        shiftmax = 0.9

        corrmasksize = Observable(1.)
        xcorr_unsm = @lift unsmoothed_xcorr_ccs(CorrelationKind.Raw(),
            mask_ccs(images_with_data[2].mod, $corrmasksize*u"mas"),
            mask_ccs(images_with_data[1].mod, $corrmasksize*u"mas");
            radecs=(@p range(0u"mas"±shiftmax*u"mas", length=64) grid(ra=reverse(__), dec=__)));
        xcorrimg = @lift smooth($xcorr_unsm; σ=$σobs*u"mas")

        curpos = Observable(SVector(NaN, NaN))
        curpos_fplt = @lift FPlot($curpos, identity; color=(:white, 0.5), strokewidth=1, strokecolor=Ref(:black), marker=marker_lw(:cross, 0.5), markersize=15)

        shiftax = Axis(fig[1,1], limits=(-1, 1, -1, 1) .* shiftmax, aspect=DataAspect(), xreversed=true, xlabel="ra (mas)", ylabel="dec (mas)", backgroundcolor=:black)
        image!(xcorrimg, colormap=:turbo, colorrange=@lift (0.6, 1) .* maximum(abs, $xcorrimg))
        text!((0.99, 0.99), text="Image shift\nColor = X-corr", space=:relative, align=(:right, :top), color=:white)
        vlines!(0, color=:white)
        hlines!(0, color=:white)
        scatter!(shiftobs; color=:white, strokecolor=:black, strokewidth=1)
        beampoly!(curbeam, color=(:white, 0), strokecolor=(:white, 0.3), strokewidth=1)

        ax_si, plt_si = image(fig[1,3], SI, colormap=:turbo, colorrange=(-1.5, 1), axis=(xlabel="", ylabel=""), interpolate=true)
        text!((0.99, 0.99), text="Spectral index", space=:relative, align=(:right, :top))
        scatter!(curpos_fplt)
        Colorbar(fig[1,4], plt_si)

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

        on(events(fig).mouseposition, priority=10) do pos
            if is_mouseinside(ax) || is_mouseinside(ax_si)
                curpos[] = Tuple(mouseposition(is_mouseinside(ax) ? ax : ax_si))
                return Consume()
            else
                curpos[] = (NaN, NaN)
            end
        end

        sg = SliderGrid(fig[0,3],
            (label = "Beam FWHM", range=0.05:0.01:1.5, format = "{:.2f} mas", startvalue=InterferometricModels.σ_to_fwhm(σobs[])),
            (label = "X-corr mask radius", range=0:0.1:3, format = "{:.2f} mas", startvalue=corrmasksize[]),
            tellheight=true, tellwidth=false)
        on(sg.sliders[1].value) do fwhm
            σobs[] = InterferometricModels.fwhm_to_σ(fwhm)
        end
        notify(sg.sliders[1].value)
        map!(identity, corrmasksize, sg.sliders[2].value)
        notify(sg.sliders[2].value)

        let hl_strs = highlight_string_differences(files; colors=[:red, :blue])
            Label(fig[0,1:2], rich(hl_strs[1] * "\n" * hl_strs[2], color=:gray70, font="JuliaMono"))
        end

        resize_to_layout!()

        MakieExtra.show_gl_icon_in_dock()
        display(fig)
    end
end


end

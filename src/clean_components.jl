module CorrelationKind
struct Raw end
struct Normalized end
end

mask_ccs(mod::MultiComponentModel, radius) =
    @modify(components(mod)) do ccs
        filter(ccs) do cc
            norm(coords(cc)) > radius
        end
    end

aggregate_ccs_by_pixel(ccs; step) =
    @modify(components(ccs)) do ccs
        @p let
            ccs
            group_vg(round_.(coords(_), step))
            map(InterferometricModels.Point(sum(flux, _), key(_)))
        end
    end

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

xcorr_ccs(kind::CorrelationKind.Raw, ccs_a, ccs_b; radecs, σ) = smooth!(unsmoothed_xcorr_ccs(kind, ccs_a, ccs_b; radecs); σ)

function cc_to_image(cc; σ, radecs, shift=SVector(0u"mas", 0u"mas"))
    img = map(Returns(0.), radecs)
    for c in components(cc)
        coo = coords(c) + shift
        img(Near(coo[1]), Near(coo[2]), :) .+= c.flux
    end
    return smooth!(img; σ)
end

spectral_index(a::MultiComponentModel, b::MultiComponentModel, freqs; shift, radecs, σ, kwargs...) = 
    spectral_index(cc_to_image(a; σ, radecs, shift), cc_to_image(b; σ, radecs), freqs; kwargs...)

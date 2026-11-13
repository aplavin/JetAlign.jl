round_(x::Number, precision::Number) = round(Integer, x / precision) * precision

"`model` without the components within `radius` of the origin."
mask_ccs(model::MultiComponentModel, radius) =
    @modify(components(model)) do ccs
        filter(cc -> norm(coords(cc)) > radius, ccs)
    end

"`model` with components merged per cell of size `step`: one point per cell with their total flux."
aggregate_ccs_by_pixel(model; step) =
    @modify(components(model)) do ccs
        @p let
            ccs
            group_vg(round_.(coords(_), step))
            map(InterferometricModels.Point(sum(flux, _), key(_)))
        end
    end

"Index `(i, j)` of the `radecs` grid pixel containing position `coo`, or `nothing` outside the grid."
function pixel_index(radecs, coo)
    ij = map(axiskeys(radecs), Tuple(coo)) do ks, x
        i = round(Int, (x - first(ks)) / step(ks)) + 1
        i ∈ eachindex(ks) ? i : nothing
    end
    any(isnothing, ij) ? nothing : ij
end

"`f(c)` of each component `c` summed into the `radecs` pixel at `coo(c)`; components outside the grid are dropped."
function grid_sum(f, coo, ccs, radecs)
    img = map(Returns(0.0), radecs)
    for c in ccs
        ij = pixel_index(radecs, coo(c))
        isnothing(ij) || (img[ij...] += f(c))
    end
    img
end

"""
    xcorr_ccs(a, b; radecs)

Cross-correlation of clean-component models `a` and `b` on the `radecs` grid of offsets, unsmoothed.
It peaks at the shift that aligns `b` with `a`.
"""
function xcorr_ccs(a, b; radecs)
    step = @p axiskeys(radecs) map(abs(Base.step(_))) uniqueonly
    pairs = Iterators.product(components.(aggregate_ccs_by_pixel.((a, b); step=step / 5))...)
    grid_sum(((ca, cb),) -> flux(ca) * flux(cb), ((ca, cb),) -> coords(ca) - coords(cb), pairs, radecs)
end

"""
    cc_to_image(model; σ, radecs, shift)

Image of clean-component `model` shifted by `shift`, on the `radecs` grid, restored with a circular
Gaussian beam of width `σ` (see [`smooth`](@ref)).
"""
cc_to_image(model; σ, radecs, shift=SVector(0, 0)u"mas") =
    smooth(grid_sum(flux, c -> coords(c) + shift, components(model), radecs); σ)

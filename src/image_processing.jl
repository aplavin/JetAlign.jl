"""
    smooth(img::KeyedArray; σ)

Convolve `img` with a circular Gaussian of width `σ`, normalized to unit peak: a point source of flux
`F` becomes a Gaussian with peak `F` (Jy/beam). Everything outside the grid counts as zero.
"""
function smooth(img::KeyedArray; σ)
    kernel(ks) = exp.(-(ks .- ks') .^ 2 ./ (2σ^2))   # separable: K_ra * img * K_dec'
    ra, dec = axiskeys(img)
    @set keyless_unname(img) = kernel(ra) * keyless_unname(img) * kernel(dec)'
end

"""
    spectral_index(a, b, freqs; dynrange)

Spectral index `α` (with `S ∝ νᵅ`) between images `a` at `freqs[1]` and `b` at `freqs[2]`.
`NaN` wherever either image is fainter than its peak divided by `dynrange`.
"""
function spectral_index(a, b, freqs; dynrange)
    amin, bmin = maximum(a) / dynrange, maximum(b) / dynrange
    lnν = log(freqs[2] / freqs[1])
    map((x, y) -> x < amin || y < bmin ? NaN : log(y / x) / lnν, a, b)
end

"""
    ellipse(a, b, pa; n=65)

Closed polygon (`n` points, the last repeating the first) of an ellipse centered at zero, with semi-axes
`a` along position angle `pa` and `b` across it. `pa` is counted from north (+dec) towards east (+ra);
points are `SVector(ra, dec)`.
"""
function ellipse(a, b, pa; n=65)
    major = SVector(sincos(pa)...)        # (sin pa, cos pa): unit vector along pa in (ra, dec)
    minor = SVector(major[2], -major[1])
    map(range(0, 2π, length=n)) do θ
        s, c = sincos(θ)
        a * c * major + b * s * minor
    end
end

"FWHM ellipse of `beam` centered at zero, see [`ellipse`](@ref)."
beampoly(beam::Beam; kwargs...) = let c = EllipticGaussian(ModelComponent(beam))
    ellipse(fwhm_max(c) / 2, fwhm_min(c) / 2, position_angle(c); kwargs...)
end

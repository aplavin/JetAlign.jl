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

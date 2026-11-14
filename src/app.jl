const Limits = @NamedTuple{ra::ClosedInterval{Float64}, dec::ClosedInterval{Float64}}

"""
All user-adjustable state: what Reset restores and what is persisted per file pair.
Angular quantities are plain numbers in mas.
"""
@kwdef struct UIState
    fwhm::Float64                           # restoring beam FWHM
    dynrange::Float64 = 1e4                 # contour floor and spectral index mask: peak / dynrange
    maskradius::Float64                     # components within this radius are excluded from the x-corr
    npix::Int = 512
    halffov::Float64
    shift::SVector{2,Float64} = SVector(0.0, 0.0)   # applied to the lower-frequency image
    sicolorrange::NTuple{2,Float64} = (-1.8, 1.2)
    shiftlimits::Limits                     # shift plot view
    skylimits::Limits                       # Stokes I and spectral index plot view
end

const NPIX_CHOICES = sort([2 .^ (5:10); 3 .* 2 .^ (4:8)])  # 32, 48, 64, 96, ..., 768, 1024

# spectral index from the single current shift, or the spread over the multishift points
struct SingleSI end
struct MultiSI end

"One image's frequency, clean-component model (fluxes in Jy), restoring beam and half field of view (mas)."
function load_image(file)
    fimg = VLBI.load(file, read_data=false)
    (;
        file,
        freq=VLBI.frequency(fimg),
        model=@modify(ustrip, VLBI.load(MultiComponentModel, file) |> RecursiveOfType(ModelComponent) |> flux),
        beam=VLBI.load(Beam, file),
        halffov=ustrip(u"mas", maximum(abs, named_axiskeys(fimg).ra)),
    )
end

skygrid(halfwidth, n) = let ks = range(-halfwidth, halfwidth, length=n) * u"mas"
    grid(ra=ks, dec=ks)
end

mutable struct App{I,S}
    const files::Vector{String}             # as given, in the order stored in shifts files and the state key
    const images::I                         # by increasing frequency
    const colors::NTuple{2,Colorant}        # per image
    const xcorr_halfwidth::Float64
    const defaults::UIState
    const store::S                          # StateStore, or nothing to not persist
    const savekey::String                   # name of the file Save shifts writes to the working directory
    ui::UIState
    shifts::Vector{SVector{2,Float64}}      # multishift points; replaced on change, never mutated
    simode::Union{SingleSI,MultiSI}
    saved::UIState                          # last persisted state
    savetime::Float64
    cursor::SVector{2,Float64}              # mouse position over the sky plots, NaN when elsewhere
    const jobs::@NamedTuple{maps::AsyncJob, xcorr::AsyncJob}
    const sicolors::CacheCell
end

function App(files; shiftsfile=nothing, colors=(colorant"red", colorant"blue"), persiststate=true)
    length(files) == 2 || error("Exactly two files must be provided for alignment, got $(length(files))")
    files = map(expanduser, files)
    all(isfile, files) || error("Files not found:\n$(join(filter(!isfile, files), '\n'))")
    images = sort(map(load_image, files), by=im -> im.freq)

    fwhm = @p images map(fwhm_average(_.beam)) geomean() ustrip(u"mas", __)
    shiftmax = 1.5fwhm
    halffov = minimum(im -> im.halffov, images)
    defaults = UIState(;
        fwhm, maskradius=1.5fwhm, halffov,
        shiftlimits=(ra=-shiftmax..shiftmax, dec=-shiftmax..shiftmax),
        skylimits=(ra=-halffov..halffov, dec=-halffov..halffov),
    )
    store = persiststate ? open_store(joinpath(@get_scratch!("cache"), "cache.db"), files) : nothing
    ui = restore(store, defaults)
    shifts = isnothing(shiftsfile) ? SVector{2,Float64}[] : load_shifts(expanduser(shiftsfile), files)
    App(files, images, Tuple(colors), 4shiftmax, defaults, store,
        string(UUIDs.uuid4())[1:8],
        ui, shifts, SingleSI(), ui, time(), SVector(NaN, NaN),
        (maps=AsyncJob(), xcorr=AsyncJob()), CacheCell())
end

reset!(app::App) = (app.ui = app.defaults; app.simode = SingleSI(); app)

set_shift!(app::App, shift) = (app.ui = setproperties(app.ui; shift); app.simode = SingleSI(); app)

set_shifts!(app::App, shifts) = (app.shifts = shifts; app.simode = isempty(shifts) ? SingleSI() : MultiSI(); app)

shifts_path(app::App) = joinpath(pwd(), "$(app.savekey).toml")

function save_shifts!(app::App)
    save_shifts(shifts_path(app), app.files, app.shifts)
    @info "Saved shifts to $(shifts_path(app))"
end

"Persist `app.ui` if it changed and the last save was at least `interval` seconds ago."
function persist!(app::App; interval)
    (app.ui == app.saved || time() - app.savetime < interval) && return
    save!(app.store, app.ui)
    app.saved = app.ui
    app.savetime = time()
end


maps_inputs(app::App) = (; app.ui.fwhm, app.ui.npix, app.ui.halffov, app.ui.shift, app.ui.dynrange, app.simode, app.shifts)

"Both Stokes I images, their contour levels, and the spectral index map."
function compute_maps(images, inputs)
    (; fwhm, npix, halffov, shift, dynrange, simode, shifts) = inputs
    σ = fwhm_to_σ(fwhm) * u"mas"
    radecs = skygrid(halffov, npix)
    imgs = (
        cc_to_image(images[1].model; σ, radecs, shift=shift * u"mas"),
        cc_to_image(images[2].model; σ, radecs),
    )
    freqs = map(im -> im.freq, images)
    si = spectral_index_map(simode, images[1].model, imgs, freqs; σ, radecs, dynrange, shifts)
    (;
        extent=(-halffov..halffov, -halffov..halffov),
        imgs=map(keyless_unname, imgs),
        peaks=map(maximum, imgs),
        levels=map(img -> maprange(log, maximum(img) / dynrange, 0.9maximum(img), length=15), imgs),
        si=keyless_unname(si),
    )
end

spectral_index_map(::SingleSI, model1, imgs, freqs; dynrange, kwargs...) = spectral_index(imgs..., freqs; dynrange)

# per pixel: the range of spectral indices over all multishift points
function spectral_index_map(::MultiSI, model1, imgs, freqs; σ, radecs, dynrange, shifts)
    sis = map(shifts) do s
        spectral_index(cc_to_image(model1; σ, radecs, shift=s * u"mas"), imgs[2], freqs; dynrange)
    end
    map((vs...) -> ClosedInterval(extrema(vs)...), sis...)
end

# per pixel: color of the mean, opacity decreasing with the spread (both clamped to the colorrange)
function si_color((lo, hi), int::ClosedInterval)
    a, b = clamp.(endpoints(int), lo, hi)
    t = ((a + b) / 2 - lo) / (hi - lo)
    isnan(t) ? RGBA{Float64}(0, 0, 0, 0) : RGBA(get(ColorSchemes.turbo, t), clamp(1 - (b - a) / 1.5, 0, 1))
end

xcorr_inputs(app::App) = (; app.ui.maskradius, app.ui.fwhm)

"Cross-correlation of the two images' clean components outside the mask radius, restored with the current beam."
function compute_xcorr(images, halfwidth, inputs)
    (; maskradius, fwhm) = inputs
    masked = map(im -> mask_ccs(im.model, maskradius * u"mas"), images)
    xc = smooth(xcorr_ccs(masked[2], masked[1]; radecs=skygrid(halfwidth, 386)); σ=fwhm_to_σ(fwhm) * u"mas")
    (; extent=(-halfwidth..halfwidth, -halfwidth..halfwidth), img=keyless_unname(xc), colorrange=(0.6, 1) .* maximum(abs, xc))
end

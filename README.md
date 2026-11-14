<p align="center">
  <img src="assets/icon.png" alt="JetAlign icon" width="160">
</p>

<h1 align="center">JetAlign.jl</h1>

<p align="center">
  An interactive tool for aligning multi-frequency jet images, aimed primarily at Very Long Baseline Interferometry (VLBI) observations.
</p>

https://github.com/user-attachments/assets/22e40d12-5ef1-4389-b9d8-3d27a028b378


## Features

- **Interactive Jet Alignment**: Semi-manual alignment of multi-frequency VLBI jet images with intuitive mouse controls
- **Spectral Index Visualization**: Compute and display spectral index maps to evaluate alignment and study jet emission properties
- **Cross-Correlation Analysis**: Image cross-correlation tools for guidance in alignment

## Installation

1. **Install Julia**, tested on Julia version 1.13

2. **Install this app**:
   ```bash
   julia -e 'using Pkg; Pkg.Apps.add("JetAlign")'
   ```
   This creates the `jetalign` command in `~/.julia/bin`; add that directory to your `PATH` if it is not there yet.

Update with `julia -e 'using Pkg; Pkg.Apps.update("JetAlign")'`, uninstall with `julia -e 'using Pkg; Pkg.Apps.rm("JetAlign")'`.

## Usage

```bash
# align two FITS images:
jetalign path/to/low_freq_image.fits path/to/high_freq_image.fits

# start from shifts saved earlier with the "Save shifts" button:
jetalign --shifts saved_shifts.toml low_freq_image.fits high_freq_image.fits
```

If you don't have FITS images and just want to explore JetAlign, you can use these example files:
```bash
wget https://astrogeo.org/images/J0433+0521/J0433+0521_S_2016_11_30_pet_map.fits
wget https://astrogeo.org/images/J0433+0521/J0433+0521_X_2016_11_30_pet_map.fits
jetalign J0433+0521_S_2016_11_30_pet_map.fits J0433+0521_X_2016_11_30_pet_map.fits
```

### Without installing the app

From a clone of the repository, the same command line runs as `julia -m JetAlign`:

```bash
git clone https://github.com/aplavin/JetAlign.jl && cd JetAlign.jl
julia --project -e 'using Pkg; Pkg.instantiate()'     # once

julia --project -t auto -m JetAlign low_freq_image.fits high_freq_image.fits
```

## Contributing

JetAlign.jl is under development. Contributions, bug reports, and feature requests are welcome!

# JetAlign.jl 🚀✨

An interactive tool for aligning multi-frequency jet images, aimed primarily at Very Long Baseline Interferometry (VLBI) observations.

https://github.com/user-attachments/assets/b8535ac4-cb6f-4344-b6f0-4144ec11eb67


## Features ✨

- **Interactive Jet Alignment** 🎯: Semi-manual alignment of multi-frequency VLBI jet images with intuitive mouse controls
- **Spectral Index Visualization** 📊: Compute and display spectral index maps to evaluate alignment and study jet emission properties
- **Cross-Correlation Analysis** 🔍: Image cross-correlation tools for guidance in alignment

## Installation & Setup 🚀

1. **Clone the Repository** 📥: Download JetAlign.jl to your preferred location
   ```bash
   git clone https://github.com/aplavin/JetAlign.jl.git
   cd JetAlign.jl
   ```

2. **Install Julia** 💎: Ensure you have Julia 1.11 installed

3. **Initialize the Environment** 🔧: Start Julia with the project environment
   ```bash
   julia --project
   ```

4. **Install Dependencies** 📦: In the Julia REPL, run:
   ```julia
   julia> using Pkg; Pkg.instantiate()
   ```

## Quick Start 🚀

Load JetAlign and start aligning your jet images:

```julia
# Load the package
julia> using JetAlign

# Basic Usage - Align two FITS images
julia> imalign_interactive([
    "path/to/low_freq_image.fits",
    "path/to/high_freq_image.fits"
])

# if you don't have FITS images and just want to explore JetAlign, use these example files:
julia> imalign_interactive([
    download("https://astrogeo.org/images/J0433+0521/J0433+0521_S_2016_11_30_pet_map.fits"),
    download("https://astrogeo.org/images/J0433+0521/J0433+0521_X_2016_11_30_pet_map.fits"),
])
```

## Contributing 🤝

JetAlign.jl is under development. Contributions, bug reports, and feature requests are welcome!

# Plasm

[![Build Status](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml?query=branch%3Amain)

Links:
- https://www.codeconvert.ai/python-to-julia-converter (good to create the first skeleton, but problems with 0-1 indices and list concatenation)


# Clone and enable dev

```bash
git clone git@github.com:scrgiorgio/Plasm.jl.git
cd Plasm.jl

julia
using Pkg
Pkg.develop(PackageSpec(path = pwd()))
exit()
```

# Run tests

```bash
julia test/viewer.jl
julia test/hpc.jl
julia test/temple.jl
julia test/manhattan.jl
julia test/properties.jl
julia test/fenvs.jl
julia test/lar.jl
julia test/arrangements2d.jl
```

# Jupyter

```
julia
using Pkg
add IJulia
add MeshCat
add CoordinateTransformations 
add Rotations 
add GeometryBasics
add Colors
add DataStructures
using IJulia
IJulia.notebook()
``````
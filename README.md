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
```

# Run tests

```bash
julia test/viewer.jl
julia test/hpc.jl
julia test/fenvs.jl
julia test/temple.jl
```


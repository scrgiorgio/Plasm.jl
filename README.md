# Plasm

[![Build Status](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml?query=branch%3Amaster)


# Run tests

```bash

# core
julia test/viewer.jl
julia test/hpc.jl
julia test/temple.jl
julia test/manhattan.jl
julia test/properties.jl
julia test/fenvs.jl

# LAR part
julia test/lar.jl
julia test/arrange2d.jl
julia test/arrange3d.jl
julia test/arrange.jl
julia test/complex.jl
```

To test notebooks:
- install Microsoft Visual Studio
- install Julia Extension
- Open `notebooks/examples.ipynb`

For Jupyter Notebooks alternatives see [ :](https://marketsplash.com/julia-ides/)


# Developing Plasm.jl

See:
- https://julialang.org/contribute/developing_package/
- activate .


```bash
julia 

# Go to the package mode
] 

# activate environment in current directory
activate .

add Combinatorics GLFW ModernGL PyCall StaticArrays Test LinearAlgebra DataStructures SparseArrays NearestNeighbors Triangulate IntervalTrees QHull CoordinateTransformations Rotations GeometryBasics Colors MeshCat FileIO MeshIO Meshing IJulia 

# exit package MODE
# CTRL+C 

using Pkg
Pkg.resolve()

exit()
```


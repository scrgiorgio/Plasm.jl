# Plasm

[![Build Status](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml?query=branch%3Amaster)

Per alberto:

```bash

# book
git checkout main

# Boolean operation
git checkout fixing-bool
```

## Run tests

```bash
julia --project=. test/hpc.jl

# lar part
julia --project=. test/bbox.jl
julia --project=. test/points.jl
julia --project=. test/plane.jl
julia --project=. test/lar.jl

# 2d arrangement
# TODO: togliere TRIANGULATE
julia --project=. test/randomlines.jl
julia --project=. test/randombubbles.jl

# 3d arrangement
julia --project=. test/twocubes.jl
julia --project=. test/randomcubes.jl
julia --project=. test/mechanical_piece.jl

# boolean
# TODO UNION e' outer atom
julia --project=. test/boolean.jl

```

## Developing Plasm.jl

Links:
- https://discourse.julialang.org/t/pattern-for-activating-the-current-project/15379/2

**Always remember to activate the current project** (in the current directory):

```bash
julia
using Pkg
Pkg.status()
Pkg.activate(".")
# Pkg.add(name="your-dependency-name-here", version="x.y.z"))
Pkg.instantiate()
Pkg.status()

exit()

# from here use
julia --project=. whatever...
```

```bash

# clone the current repository
cd ~
mkdir -p github.com/scrgiorgio
cd github.com/scrgiorgio
git clone https://github.com/scrgiorgio/Plasm.jl
cd Plasm.jl

julia

# Go to the package mode
# see https://julialang.org/contribute/developing_package
] 

# Activate the environment in the current directory
activate .

# add packages
add Combinatorics GLFW ModernGL PyCall StaticArrays Test LinearAlgebra DataStructures SparseArrays NearestNeighbors Triangulate IntervalTrees CoordinateTransformations Rotations GeometryBasics Colors MeshCat FileIO MeshIO Meshing IJulia 

# exit package MODE
# CTRL+C 

using Pkg
Pkg.resolve()
```

## (OPTIONAL) Julia notebooks

To test notebooks:
- install Microsoft Visual Studio
- install Julia Extension
- Open `notebooks/examples.ipynb`

For Jupyter Notebooks alternatives see [ :](https://marketsplash.com/julia-ides/)

```bash
using IJulia
notebook()
# jupyterlab() (OPTIONAL) if you want to install lab

# force to use julia internal jupyter
ENV["JUPYTER"]=""
Pkg.build("IJulia")

exit()

# for Alberto's laptop: use `lab` instead of notebook
# ~/.julia/conda/3/bin/jupyter notebook --ip='*' --NotebookApp.token='' --NotebookApp.password=''
```


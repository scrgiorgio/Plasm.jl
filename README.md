# Plasm

[![Build Status](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml?query=branch%3Amaster)


This repo:
- https://github.com/scrgiorgio/Plasm.jl



## Run tests

```bash
# === geometry ===
julia --project=.  test/bbox.jl
julia --project=.  test/points.jl
julia --project=.  test/plane.jl

# === hpc === 
julia --project=.  test/hpc.jl
julia --project=.  test/structural_frames.jl

# === lar === 
julia --project=.  test/lar.jl
julia --project=.  test/random-tetgen.jl
julia --project=.  test/simplexn.jl
julia --project=.  test/arrange2d.jl
julia --project=.  test/arrange3d.jl
julia --project=.  test/boolean.jl
julia --project=.  test/building.jl
julia --project=.  test/svg.jl

```


## Updating to Latest Version

After `git pull`:

```bash
# Update dependencies
julia --project=. -e "using Pkg; Pkg.resolve(); Pkg.instantiate()"

# Test the installation
julia --project=. -e "using Plasm; VIEW(CUBE(1))"
```

## Developing Plasm.jl

Links:
- https://discourse.julialang.org/t/pattern-for-activating-the-current-project/15379/2

**Always remember to activate the current project** (in the current directory):

```bash

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

using Pkg
Pkg.activate(".")

Pkg.add([
  "Combinatorics", 
  "CSV",
  "DataFrames",
  "GLFW", 
  "ModernGL", 
  "PyCall", 
  "StaticArrays", 
  "DataStructures", 
  "NearestNeighbors", 
  "Triangulate", 
  "IntervalTrees",
  "CoordinateTransformations", 
  "Rotations",
  "GeometryBasics",
  "Colors",
  "FileIO",
  "MeshCat",
  "IJulia",
  "PlyIO",
  "XML",
  "JSON"
  ])

# Standard libraries (no need to add, already included with Julia):
# using Test, LinearAlgebra, Random, Statistics, SparseArrays

# update the manifest too
Pkg.resolve()

exit()
```

## Viewer options

- GLMakie.jl (using internally glfw for desktop)
  - more complicate to use
  - more development by Julia Communinity
  - HUGE problem: takes forever to run. Not a viable option (!)
    - see https://discourse.julialang.org/t/starting-glmakie-takes-very-long/64106/2

- MeshCat.jl (using MeshCat/ThreeJL)
  - no support for 3d text (we can always use alberto's code)
  - no many releases
  - seems the best choice

## (OPTIONAL) Julia notebooks

To test notebooks:
- install Microsoft Visual Studio
- install Julia Extension
- Open `notebooks/examples.ipynb`

For Jupyter Notebooks alternatives see [ :](https://marketsplash.com/julia-ides/)

```bash
using IJulia
notebook(dir=pwd())
# jupyterlab() (OPTIONAL) if you want to install lab

# force to use julia internal jupyter
ENV["JUPYTER"]=""
Pkg.build("IJulia")

exit()

# for Alberto's laptop: use `lab` instead of notebook
# ~/.julia/conda/3/bin/jupyter notebook --ip='*' --NotebookApp.token='' --NotebookApp.password=''
```


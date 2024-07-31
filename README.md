# Plasm

[![Build Status](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml?query=branch%3Amaster)


## Run tests

```bash

# core
julia --project=. test/viewer.jl
julia --project=. test/hpc.jl
julia --project=. test/temple.jl
julia --project=. test/manhattan.jl
julia --project=. test/properties.jl
julia --project=. test/fenvs.jl

# LAR part
julia --project=. test/lar.jl
julia --project=. test/complex.jl

# deprecated
# julia --project=. test/arrange2d.jl
# julia --project=. test/arrange3d.jl
# julia --project=. test/arrange.jl

julia --project=. test/randomlines.jl
julia --project=. test/randomcubes.jl


julia --project=. test/congruence.jl
julia --project=. test/3tubes3d.jl
```


## Developing Plasm.jl

Links:
- https://discourse.julialang.org/t/pattern-for-activating-the-current-project/15379/2

**Always remember to activate the current project** (in the current directory):

```bash
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
add Combinatorics GLFW ModernGL PyCall StaticArrays Test LinearAlgebra DataStructures SparseArrays NearestNeighbors Triangulate IntervalTrees QHull CoordinateTransformations Rotations GeometryBasics Colors MeshCat FileIO MeshIO Meshing IJulia 

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


# Plasm

[![Build Status](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml?query=branch%3Amaster)


# Run tests

```bash

# alberto's laptop
# alias julia=julia19

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

```bash

# clone the current repository
cd ~
mkdir -p github.com/scrgiorgio
cd github.com/scrgiorgio
git clone https://github.com/scrgiorgio/Plasm.jl
cd Plasm.jl

# alberto's laptop
# alias julia=julia19

julia 

# Go to the package mode
# see https://julialang.org/contribute/developing_package
] 

# Activate the environment in the current directory
activate .

add Combinatorics GLFW ModernGL PyCall StaticArrays Test LinearAlgebra DataStructures SparseArrays NearestNeighbors Triangulate IntervalTrees QHull CoordinateTransformations Rotations GeometryBasics Colors MeshCat FileIO MeshIO Meshing IJulia 

# exit package MODE
# CTRL+C 

using Pkg
Pkg.resolve()

# install jupyter notebook
using IJulia
notebook()
# jupyterlab() (OPTIONAL) if you want lab

# force to use julia internal jupyter
ENV["JUPYTER"]=""
Pkg.build("IJulia")

exit()

git diff

# for Alberto's laptop: use `lab` instead of notebook
~/.julia/conda/3/bin/jupyter notebook --ip='*' --NotebookApp.token='' --NotebookApp.password=''
```


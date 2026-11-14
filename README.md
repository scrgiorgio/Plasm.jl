# Plasm

[![Build Status](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml?query=branch%3Amain)


This repo:
- https://github.com/scrgiorgio/Plasm.jl


## Setup

Requires Julia >= 1.10.  From the repository root:

```bash
julia --version                                                        # >= 1.10
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"  # no ✗ errors
julia --project=. -e "using Plasm; VIEW(CUBE(1))"                      # opens the viewer
```


## Tutorials

```bash
julia --project=. tutorials/point-cloud.jl ./tutorials/autzen_trim.laz 
julia --project=. tutorials/point-cloud.jl ./Roma_LIDAR.laz --color-by-height
```


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


## Updating

After `git pull` (new packages):

```bash
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
```

After a Julia upgrade, or on "resolved with a different julia version" / `✗` errors:

```bash
rm Manifest.toml
julia --project=. -e "using Pkg; Pkg.update(); Pkg.precompile()"
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


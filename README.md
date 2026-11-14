# Plasm

[![Build Status](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml?query=branch%3Amain)


This repo:
- https://github.com/scrgiorgio/Plasm.jl


## Setup

Requires Julia >= 1.10.  From the repository root:

```bash

# make sure julia >= 1.10
julia --version                                                        # >= 1.10

# Julia 1.10 only (BPA.jl is not registered)
julia --project=. -e "using Pkg; Pkg.add(url=\"https://github.com/ctsilva/BPA.jl\")"

# check no ✗ errors
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"

# opens the viewer
julia --project=. -e "using Plasm; VIEW(CUBE(1))" 
```


## Tutorials

### Point cloud:

Read a LAZ file and render it as a point cloud.

Options: 
- `--max-points=N` subsamples big files (default 2000000);
- `--color-by-height` colors the points by elevation (blue=low, red=high)

```bash
julia --project=. tutorials/point-cloud.jl ./tutorials/autzen_trim.laz
julia --project=. tutorials/point-cloud.jl ./Roma_LIDAR.laz --max-points=500000 --color-by-height
```

### Ball pivoting

Surface reconstruction of a LIDAR point cloud with the Ball Pivoting Algorithm ([BPA.jl](https://github.com/ctsilva/BPA.jl)).

Options:
- `--max-points=N` subsamples big files (default 200000)
- `--radii=a,b,...` ball radii as multiples of the point spacing (typical distance between neighbouring
  points, printed by the tutorial), one reconstruction pass per radius, in increasing order (default 2,4)

How to choose the radii (the tutorial prints the point spacing and tips for the next run):
- the ball must be larger than the spacing, or it falls through the points: start at 1.5-2
- it must be smaller than the details to keep, or it bridges over them (e.g. tree crowns, walls)
- for unevenly sampled data add larger radii: each pass fills holes the previous one left
- look at "points used": if many points are not used, increase the radii
- e.g. on `autzen_trim.laz`: `--radii=1.5` uses 91% of the points, `--radii=2` 93%, `--radii=2,4` 95%

In the viewer press `L` to show the triangle edges.

```bash
julia --project=. tutorials/bpa-example.jl ./tutorials/autzen_trim.laz
julia --project=. tutorials/bpa-example.jl ./Roma_LIDAR.laz --max-points=300000 --radii=2,4
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


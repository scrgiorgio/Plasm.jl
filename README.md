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

Read a LAZ file and render it as a point cloud, optionally stored in an octree.

Options: 
- `--max-points=N` subsamples big files (default 2000000);
- `--color-by-height` colors the points by elevation (blue=low, red=high)
- `--octree-depth=d` stores the points in an octree with 8^d leaves, each leaf with its own color (default 0, no octree)
- `--octree-box=xmin,xmax,ymin,ymax,zmin,zmax` octree root box, min and max for each axis, in the file coordinates
  (default: the point cloud bounding box, printed by the tutorial in the same format); points outside are dropped

`POINTCLOUD(V, octree_depth=d, octree_box=(min, max))` stores the points in an octree:
- the root box is `octree_box`, a 3D box given by its min and max corners, e.g. `octree_box=([0.0,0,0], [100.0,50,20])`,
  or a `Box3d`/`BoxNd`; if not given, the point cloud bounding box (points outside `octree_box` are dropped with a warning)
- every node is an `Hpc` with 8 children, X varies fastest: child `k = ix + 2*iy + 4*iz + 1`
- leaves (at depth `d`) hold the points, each point in exactly one leaf; empty leaves are kept
- every node has properties `"octree_box" => (min, max)` (its 3D box corners) and `"octree_depth" => level` (root = 0)

```bash
# no octree
julia --project=. tutorials/point-cloud.jl ./tutorials/autzen_trim.laz

julia --project=. tutorials/point-cloud.jl ./Roma_LIDAR.laz --max-points=500000 --color-by-height

# octree enabled
julia --project=. tutorials/point-cloud.jl ./tutorials/autzen_trim.laz --octree-depth=3

julia --project=. tutorials/point-cloud.jl ./tutorials/autzen_trim.laz --octree-depth=2 --octree-box=636000,636600,848900,849500,400,530

julia --project=. tutorials/point-cloud.jl ./Roma_LIDAR.laz --max-points=500000 --octree-depth=4
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

# === robustness === 
julia --project=.  test/many-cubes.jl

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


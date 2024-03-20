# Plasm

[![Build Status](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml?query=branch%3Amain)


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

# unused
# julia test/organizer.jl
# import Pkg; Pkg.add("ViewerGL")` to install the ViewerGL package.
```

# View properties

Important: use String as Dict keys, not symbols

```julia
View(batches, properties=Dict(
  "background_color" => Point4d(1,1,1,1))
  "title"            => "my title")
  "use_ortho"        => false)
  "show_lines"       => true)
  "fov"              => 60.0)
  "pos"              => Point3d(0,0,0) )
  "dir"              => Point3d(0,0,-1))
  "vup"              => Point3d(1,0,0))
  "znear"            => 0.00001) 
  "zfar "            => 100.0) 
  "walk_speed "      => 10.0) 
)
'''

# Jupyter or Juno?

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
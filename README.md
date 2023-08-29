# Plasm

[![Build Status](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml?query=branch%3Amain)

Links:
- https://github.com/scrgiorgio/Plasm.jl
- https://julialang.org/contribute/developing_package/

Translators:
- https://github.com/GunnarFarneback/CrudePythonTranslator.jl  (NOT GOOD)
- https://www.codeconvert.ai/python-to-julia-converter


# Setup

```
set PATH=%PATH%;c:/Julia-1.9.3\bin
julia

# enable dev of current directory
using Pkg
Pkg.develop(PackageSpec(path = pwd()))
exit()

# run tests
julia test/viewer.jl
julia test/hpc.jl
julia test/fenvs.jl
julia test/temple.jl
```


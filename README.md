# Plasm

[![Build Status](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/scrgiorgio/Plasm.jl/actions/workflows/CI.yml?query=branch%3Amain)

Links:
- https://github.com/scrgiorgio/Plasm.jl
- https://julialang.org/contribute/developing_package/

Translators:
- https://github.com/GunnarFarneback/CrudePythonTranslator.jl  (NOT GOOD)
- https://www.codeconvert.ai/python-to-julia-converter


# Setup

```bash
set PATH=%PATH%;c:\Julia-1.9.2\bin
julia

using Pkg
Pkg.add("Combinatorics)
Pkg.add("ModernGL")
Pkg.add("GLFW")
Pkg.add("StaticArrays")
Pkg.add("PyCall")
exit()


Tests:

```bash
julia src/viewer.jl
julia src/hpc.jl
julia src/fenvs.jl
julia src/temple.jl
- 
- 
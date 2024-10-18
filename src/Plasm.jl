module Plasm

using LinearAlgebra
using DataStructures
using SparseArrays
using StaticArrays
using Triangulate
using NearestNeighbors
using Random
using Statistics

TO_GEOMETRY_DEFAULT_PRECISION_DIGITS = 14
export TO_GEOMETRY_DEFAULT_PRECISION_DIGITS

DEFAULT_LAR_FONT_SIZE=0.01
export DEFAULT_LAR_FONT_SIZE

DEFAULT_VIEWER="glfw"
export DEFAULT_VIEWER

include("./config.jl")
include("./utils.jl")
include("./geometry.jl")
include("./points.jl")
include("./plane.jl")
include("./defaults.jl")
include("./viewer.jl")
include("./hpc.jl")
include("./fenvs.jl")

LAR_DEFAULT_ERR=1e-8
export LAR_DEFAULT_ERR

LAR_ARRANGE_VERSION=1
export LAR_ARRANGE_VERSION

# make them match 
LAR_EXPERIMENTAL_ARRANGE_ROUND            =    4
LAR_EXPERIMENTAL_ARRANGE_PERTURBATION     = 1e-4 * 0.1 
LAR_ARRANGE3D_USE_EXPERIMENTAL = true

include("./lar/lar.jl")
include("./lar/simplexn.jl")
include("./lar/integr.jl")
include("./lar/classify.jl")
include("./lar/arrange2d.jl")
include("./lar/arrange3d.jl")
include("./lar/experimental/arrange2d.jl")
include("./lar/experimental/arrange3d.jl")
include("./lar/boolean.jl")


function __init__()
	InitPythonHullCode()
end

end

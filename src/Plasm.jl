module Plasm

TO_GEOMETRY_DEFAULT_PRECISION_DIGITS = 14
export TO_GEOMETRY_DEFAULT_PRECISION_DIGITS

DEFAULT_LAR_FONT_SIZE=0.04
export DEFAULT_LAR_FONT_SIZE

DEFAULT_VIEWER="glfw"
export DEFAULT_VIEWER

# scrgiorgio: do not use different errors/tolerance in lar code, try to use same number
LAR_DEFAULT_ERR=1e-8
export LAR_DEFAULT_ERR


include("./config.jl")
include("./geometry.jl")
include("./defaults.jl")

include("./viewer.jl")

include("./hpc.jl")
include("./fenvs.jl")

using LinearAlgebra
using DataStructures
using SparseArrays
using StaticArrays
using Triangulate
using NearestNeighbors
using Random
using Statistics

include("./lar/points.jl")
include("./lar/plane.jl")
include("./lar/bbox.jl")
include("./lar/lar.jl")
include("./lar/classify.jl")

# from Alberto: do not touch. Too complicate to do any reorganization right now
include("./lar/arrange2d.jl")
include("./lar/arrange3d.jl")

# scrgiorgio: to fix
include("./lar/boolean.jl")

function __init__()
	InitPythonHullCode()
end

end

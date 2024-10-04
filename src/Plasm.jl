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

DEFAULT_LAR_FONT_SIZE=0.03
export DEFAULT_LAR_FONT_SIZE

DEFAULT_VIEWER="glfw"
export DEFAULT_VIEWER

# scrgiorgio: do not use different errors/tolerance in lar code, try to use same number
LAR_DEFAULT_ERR=1e-8
export LAR_DEFAULT_ERR

LAR_FRAGMENT_ERR=1e-5
LAR_FRAGMENT_DIGITS=4
export LAR_FRAGMENT_ERR
export LAR_FRAGMENT_DIGITS

include("./config.jl")
include("./geometry.jl")
include("./points.jl")
include("./plane.jl")
include("./defaults.jl")
include("./viewer.jl")
include("./hpc.jl")
include("./fenvs.jl")
include("./lar/lar.jl")
include("./lar/classify.jl")
include("./lar/arrange2d.jl")
include("./lar/arrange3d.jl")
include("./lar/experimental/arrange2d.jl")
include("./lar/experimental/arrange3d.jl")
include("./lar/boolean.jl")

LAR_ARRANGE_VERSION=1

function ARRANGE2D(lar::Lar; debug_mode=false)
	return LAR_ARRANGE_VERSION==2 ? arrange2d_v2(lar,debug_mode=debug_mode) : arrange2d_v1(lar)
end

function ARRANGE3D(lar::Lar)
	return LAR_ARRANGE_VERSION==2 ? arrange3d_v2(lar,debug_mode=debug_mode) : arrange3d_v1(lar)
end

function INNERS(lar)
	return LAR_ARRANGE_VERSION==2 ? arrange3d_v2_inners(lar,debug_mode=debug_mode) : arrange3d_v1_inners(lar)
end

function OUTERS(lar)
	return LAR_ARRANGE_VERSION==2 ? arrange3d_v2_outers(lar,debug_mode=debug_mode) : arrange3d_v1_outers(lar)
end

export ARRANGE2D
export ARRANGE3D
export INNERS
export OUTERS


function __init__()
	InitPythonHullCode()
end

end

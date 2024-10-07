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

DEFAULT_LAR_FONT_SIZE=0.05
export DEFAULT_LAR_FONT_SIZE

DEFAULT_VIEWER="glfw"
export DEFAULT_VIEWER

LAR_DEFAULT_ERR=1e-8
export LAR_DEFAULT_ERR

include("./config.jl")
include("./utils.jl")
include("./geometry.jl")
include("./points.jl")
include("./plane.jl")
include("./defaults.jl")
include("./viewer.jl")
include("./hpc.jl")
include("./fenvs.jl")
include("./lar/lar.jl")
include("./lar/classify.jl")

LAR_ARRANGE2D_ROUND         = 4
LAR_ARRANGE2D_PERTURBATION  = 1e-4 * 0.01

LAR_ARRANGE2D_SMALL_TRIANGLES_ERR    = 1e-4
LAR_ARRANGE3D_UNPROJECT_ROUND_DIGITS = 4

include("./lar/arrange2d.jl")
include("./lar/arrange3d.jl")


include("./lar/experimental/arrange2d.jl")
include("./lar/experimental/arrange3d.jl")

include("./lar/boolean.jl")

LAR_ARRANGE_VERSION=1

function ARRANGE2D(lar::Lar; debug_mode=false)::Lar
	return LAR_ARRANGE_VERSION==2 ? arrange2d_v2(lar) : arrange2d_v1(lar)
end

function ARRANGE3D(lar::Lar; debug_mode=false)::Lar
	return LAR_ARRANGE_VERSION==2 ? arrange3d_v2(lar,debug_mode=debug_mode) : arrange3d_v1(lar)
end

function INNERS(lar::Lar)::Lar
	return LAR_ARRANGE_VERSION==2 ? arrange3d_v2_inners(lar) : arrange3d_v1_inners(lar)
end

function OUTERS(lar::Lar)::Lar
	return LAR_ARRANGE_VERSION==2 ? arrange3d_v2_outers(lar) : arrange3d_v1_outers(lar)
end

export LAR_ARRANGE_VERSION
export ARRANGE2D
export ARRANGE3D
export INNERS
export OUTERS

function __init__()
	InitPythonHullCode()
end

end

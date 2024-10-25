module Plasm

using LinearAlgebra
using DataStructures
using SparseArrays
using StaticArrays
using Triangulate
using NearestNeighbors
using Random
using Statistics

include("./config.jl")
include("./utils.jl")
include("./geometry.jl")
include("./points.jl")
include("./plane.jl")
include("./defaults.jl")

TEXT_LINE_WIDTH=2
DEFAULT_VIEWER="glfw"
include("./viewer.jl")

TO_GEOMETRY_PRECISION = 1e-14
include("./hpc.jl")

MAP_PRECISION=1e-8
include("./fenvs.jl")

# ///////////////////////////////////////////////

LAR_DEFAULT_ERR=1e-8

LAR_VERTEX_FONT_SIZE = 0.05
LAR_EDGE_FONT_SIZE   = 0.06
LAR_FACE_FONT_SIZE   = 0.07

LAR_VERTEX_COLOR     = BLACK
LAR_EDGE_COLOR       = BLUE
LAR_FACE_COLOR       = MAGENTA

# not used right now
LAR_EXPERIMENTAL_ARRANGE_PERTURBATION = 1e-4 * 0.1 

# use new TGW and SPLIT code (cannot be sure it's more sbable, need testing)
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

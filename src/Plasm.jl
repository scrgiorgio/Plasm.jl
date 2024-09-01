module Plasm

include("./config.jl")
include("./viewer.jl")
include("./hpc.jl")
include("./fenvs.jl")

using LinearAlgebra
using DataStructures
using SparseArrays
using StaticArrays
using Triangulate
using NearestNeighbors

# ///////////////////////////////////////////////////////////
const LAR_DEFAULT_ERR=1e-8

include("./lar/points.jl")      # OK
include("./lar/bbox.jl")        # OK
include("./lar/dense.jl")       # OK
include("./lar/sparse.jl")      # OK
include("./lar/view.jl")        # OK

# TODO
include("./lar/classify.jl")
include("./lar/triangulate.jl")
include("./lar/arrange2d.jl")
include("./lar/arrange3d.jl")
include("./lar/boolean.jl")

function __init__()
	InitPythonHullCode()
end

end

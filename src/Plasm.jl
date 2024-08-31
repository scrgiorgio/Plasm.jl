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

include("./lar/points.jl")      # OK
include("./lar/bbox.jl")        # OK
include("./lar/dense.jl")       # OK
include("./lar/sparse.jl")

include("./lar/triangulate.jl")
include("./lar/view.jl")

# TODO
include("./lar/arrange2d.jl")
include("./lar/arrange3d.jl")
include("./lar/boolean.jl")

function __init__()
	InitPythonHullCode()
end

end

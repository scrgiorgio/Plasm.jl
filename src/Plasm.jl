module Plasm

include("./config.jl")
include("./viewer.jl")
include("./hpc.jl")
include("./fenvs.jl")

using LinearAlgebra
using DataStructures
using SparseArrays
using StaticArrays

include("./lar/points.jl")
include("./lar/bbox.jl")
include("./lar/dense.jl")
include("./lar/sparse.jl")
include("./lar/triangulate.jl")
include("./lar/view.jl")

# TODO
include("./lar/arrange.jl")
include("./lar/boolean.jl")

function __init__()
	InitPythonHullCode()
end

end

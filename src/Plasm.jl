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

# scrgiorgio: do not use different errors/tolerance in lar code, try to use same number
const LAR_DEFAULT_ERR=1e-8
export LAR_DEFAULT_ERR

include("./lar/points.jl")      
include("./lar/bbox.jl")        
include("./lar/dense.jl")       
include("./lar/sparse.jl")      
include("./lar/view.jl")        
include("./lar/classify.jl")    

# from Alberto: do not touch. Too complicate to do any reorganization right now
include("./lar/triangulate.jl")
include("./lar/arrange2d.jl")
include("./lar/arrange3d.jl")

# scrgiorgio: to fix
include("./lar/boolean.jl")

function __init__()
	InitPythonHullCode()
end

end

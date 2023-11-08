module Plasm
include("./viewer.jl")
include("./hpc.jl")
include("./fenvs.jl")

# module init
function __init__()
	InitToLAR()
end

end
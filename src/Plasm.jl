module Plasm

include("./config.jl")
include("./viewer.jl")
include("./hpc.jl")
include("./fenvs.jl")
include("./arrange.jl")
include("./complex.jl")
include("./boolean.jl")
include("./congruence.jl")
include("./lar2triangles.jl")

# module init
function __init__()
	InitToLAR()
end

end

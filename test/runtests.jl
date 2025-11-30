using Test
using Plasm

@testset "Plasm.jl" begin
    # Include all test files
    include("bbox.jl")
    #include("points.jl")
    #include("plane.jl")
    #include("hpc.jl")
    #include("structural_frames.jl")
    #include("lar.jl")
    #include("simplexn.jl")
    #include("arrange2d.jl")
    #include("arrange3d.jl")
    #include("boolean.jl")
    #include("svg_parser.jl")
end

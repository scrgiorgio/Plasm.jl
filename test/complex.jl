using Test, Plasm, LinearAlgebra

import Base.+  
+(f::Function, g::Function) = (x...) -> f(x...) + g(x...)  

# //////////////////////////////////////////////////////////////////////////////
@testset "IsSimplex" begin
   @test IsPolytope(CUBE(100))
   @test IsSimplex(CUBOID([1]))  == true
   @test IsSimplex(SIMPLEX(1))  == true
   @test IsSimplex(CUBOID([1,1])) == false
   @test IsPolytope(CUBOID([1,1]))  == true
   @test IsPolytope(CUBOID([3,3]))  == true
   @test IsPolytope(CUBE(100))
end

# //////////////////////////////////////////////////////////////////////////////
@testset "SIMPLEX" begin
   @test SIMPLEX(1).T == MatrixNd(2)
   @test SIMPLEX(4).T == MatrixNd(5)
   @test SIMPLEX(20).childs[1].hulls == [1:21]
   @test SIMPLEX(20).childs[1].points[1] == zeros(20)
end

simplexfacets

# //////////////////////////////////////////////////////////////////////////////
@testset "CUBOID" begin
   @test typeof(CUBOID([1]).childs[1].childs[1].points) == Vector{Vector{Float64}}
   @test typeof(CUBOID([1]).childs[1].childs[1].hulls) == Vector{Vector{Int64}}
   @test CUBOID([1]).childs[1].childs[1].points == [[0.0],[1.0]]
   @test CUBOID([1]).childs[1].childs[1].hulls == [[1,2]]
   @test typeof(CUBOID([1,1,1,1]).childs[1].childs[1].points[1]) == Vector{Float64}
   @test typeof(CUBOID([1,1,1,1]).childs[1].childs[1].hulls[1]) == Vector{Int64}
   @test CUBOID(DIESIS(4)(1)).childs[1].childs[1].points[1] == [0.,0.,0.,0.]
   @test CUBOID(DIESIS(4)(1)).childs[1].childs[1].hulls[1] == 1:2^4
end;
@testset "CUBOID" "$shape" for shape in [[1],[3,2,1],[3,2],[1,1,1,1]]
   @test CUBOID(shape).childs[1].childs[1].points[1] == zeros(LEN(shape))
end;

# //////////////////////////////////////////////////////////////////////////////
using QHull # REMARK

@testset "CHULL" begin
   points = rand(6, 3)
   obj = CHULL(points)
   VIEW(Hpc(obj.V, obj.C[:EV]))
end

# //////////////////////////////////////////////////////////////////////////////













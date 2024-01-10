using Test, Plasm, LinearAlgebra
using QHull 


# extracted from alberto's code

"""
SIMPLEX(1)
SIMPLEX(2)
SIMPLEX(3)
SIMPLEX(4)

CUBOID([1,1])
CUBOID([1,1,1])
CUBOID([1,1,1,1])

simplex6 = simplex(6,complex=true)
simplex6.V
simplex6.C
Hpc(simplex6.V, simplex6.C[:C5V])

Simplex6 = Lar(Hpc(simplex6.V, simplex6.C[:C5V]))

obj3 = Hpc(Lar(CUBOID([1,1,1]))); # OK
VIEW(Hpc(Lar(obj3).V, Lar(obj3).C[:EV])) # OK
VIEW(Hpc(Lar(obj3).V, Lar(obj3).C[:FV])) # OK

LAR(STRUCT(SQUARE(1), T(1,2)(0.5,0.25), SQUARE(1))) # 2D
LAR(STRUCT(CUBE(1), T(1,2)(0.5,0.25), CUBE(1))) # 3D |V|=16,|FV|=12,|EV|=24 ... OK!
LAR(STRUCT(CUBE(1), T(1)(1), CUBE(1))) # |V|=12,|FV|=11,|EV|=20 ... OK!  correct by using LAR ... uncorrect using Lar !!!

cube3 = (LAR∘STRUCT∘CAT∘DIESIS(3))([T(1)(1), CUBE(1)]); # NB:  LAR  !!!
VIEW(Hpc(cube3.V, cube3.C[:EV])) # OK
VIEW(Hpc(cube3.V, cube3.C[:FV])) # OK

grid1D = INTERVALS(10.0)(10) # equivalent to the next
grid1D = QUOTE(DIESIS(10)(1.0)) # equivalent to the former

mesh = CUBOIDGRID([5,5,5])
VIEW(Hpc(mesh.V, mesh.C[:EV]))

GRID1(5)

mesh = CUBOIDGRID([10,10])
VIEW(Hpc(mesh.V, mesh.C[:EV])) # problema ai bordi
VIEW(Hpc(mesh.V, mesh.C[:FV])) # OK ai bordi
"""


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

# simplexfacets

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
# scrgiorgio: this fails for me
# @testset "CHULL" begin
#    points = rand(6, 3)
#    obj = CHULL(points)
#    VIEW(Hpc(obj.V, obj.C[:EV]))
# end














using Test
using Plasm, LinearAlgebra



@testset "dense.jl" begin

   
   @testset "struct Lar" begin
      obj = SIMPLEX(3)
      larobj = LAR(obj)
      @test typeof(larobj) == Lar
      @test typeof(larobj.V) == Matrix{Float64}
      geo = ToGeometry(obj, precision=9)
	   n = length(geo.points)    # number of vertices  (cols of V)
	   m = length(geo.points[1]) # embedding dime (rows of V) n. of coords
      @test size(larobj.V,2) == 4
      @test size(larobj.V,1) == 3
      @test larobj.V == [
        0.0  1.0  0.0  0.0
        1.0  0.0  0.0  0.0
        0.0  0.0  0.0  1.0]
      @test larobj.C == Dict{Symbol, AbstractArray}(:CV => [[1, 2, 3, 4]], :FV => [[1, 2, 3], [2, 3, 4], [1, 3, 4], [1, 2, 4]], :EV => [[2, 3], [1, 3], [1, 2], [3, 4], [2, 4], [1, 4]])
   end;

   @testset "LAR" begin
      obj = Q(1)*Q(1)*Q(1)
      @test typeof(obj) == Hpc
      @test typeof(LAR(obj)) == Lar
      @test (LAR(obj)).V == [1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0; 
      1.0 1.0 0.0 0.0 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0]
      @test (LAR(obj)).C[:CV] == [[1, 2, 3, 4, 5, 6, 7, 8]]
      @test (LAR(obj)).C[:FV] == [[1, 2, 3, 4], [3, 4, 5, 6], [1, 3, 5, 7], 
      [2, 4, 6, 8], [1, 2, 7, 8], [5, 6, 7, 8]]
      @test (LAR(obj)).C[:EV] == [[3, 4], [2, 4], [1, 2], [1, 3], [5, 6], 
      [4, 6], [3, 5], [5, 7], [1, 7], [6, 8], [2, 8], [7, 8]]
   end;
   
   @testset "HPC" begin
      obj = Q(1)*Q(1)*Q(1)
      points = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [0.0, 1.0,
       1.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0]]
      hulls = [[1, 2, 3, 4, 5, 6, 7, 8]]
      obj = MkPol(points, hulls).childs[1]
      @test dim(obj) == 3
      @test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
      obj = ToSimplicialForm(obj)
      @test dim(obj) == 3
      @test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
      @test length(obj.points) == 8
      @test length(obj.hulls) == 6
   end;
      
   @testset "CUBOIDGRID" begin
      grid211 = CUBOIDGRID([2,1,1])
      @test typeof(grid211) == Lar
      @test size(grid211.V) == (3, 12)
      @test length(grid211.V) == 36
      @test length(grid211.C[:CV]) == 2
      @test length(grid211.C[:FV]) == 11
      @test length(grid211.C[:EV]) == 20
      grid213 = CUBOIDGRID([2,1,3])
      @test typeof(grid211) == Lar
      @test size(grid213.V) == (3, 24)
      @test length(grid213.V) == 72
      @test length(grid213.C[:CV]) == 6
      @test length(grid213.C[:FV]) == 29
      @test length(grid213.C[:EV]) == 46
   end;

   @testset "misc" begin
      hpc = Plasm.CUBE(1) 
      lar = LAR(hpc)
      V, EV, FV  = lar.V, lar.C[:EV], lar.C[:FV]
      
      # each vertex has 3 edges
      VE=compute_VE(EV)
      @test(length(VE)==8 && all([length(it)==3 for it in VE]))

      # each face has 4 edges
      FE=compute_FE(EV,FV,double_check=true)
      @test(length(FE)==6 && all([length(it)==4 for it in FE]))
   end
   

   @testset "simplify" begin
      @test(SIMPLIFY([3,2,1,1,2,3])==[1,2,3])

      lar=SIMPLIFY(Lar(BYCOL([0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0 ]),  Dict(
      :EV => [[2,1],[3,2],[1,2],[1,4],[4,3],[2,3],[3,4],[4,1]],
      :FV => [[1,2,3],[3,2,1],[4,3,1],[1,3,4]],
      :CF => [[2,1],[1,2,3,4],[1,2,3],[1,3,4],[4,3]]
      )))
      @test(lar.C[:EV]==[[1,2],[2,3],[1,4],[3,4]])
      @test(lar.C[:FV]==[[1,2,3],[1,3,4]])
      @test(lar.C[:CF]==[[1],[1,2],[1,2],[1,2],[2]])
end


end

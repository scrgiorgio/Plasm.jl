using Test, Plasm, LinearAlgebra

@testset "simplexn" begin

   @testset "LARSIMPLEX" begin
      @test LARSIMPLEX(1) == ([0.0 1.0], [[1, 2]])
      @test LARSIMPLEX(3) == ([0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0], [[1, 2, 3, 4]])
      @test length(LARSIMPLEX(5)[2]) == 1
      @test length(LARSIMPLEX(5)[2][1])==6
   end

   @testset "SIMPLEXFACETS" begin
      @test SIMPLEXFACETS([[1,2,3], [2,4,3], [1,3,5]])==[[1, 2], [1, 3], [1, 5], [2, 3], [2, 4], [3, 5], [4, 3]]
      @test SIMPLEXFACETS([[1, 2], [1, 3], [1, 5], [2, 3], [2, 4], [3, 5], [4, 3]])==[[1], [2], [3], [4], [5]]
		@test SIMPLEXGRID([1,1])[1]==[0 1 0 1; 0 0 1 1]
		@test SIMPLEXGRID([1,1])[2]==Array{Int64,1}[[1, 2, 3], [2, 3, 4]]
   end
   
	@testset "EXTRUDESIMPLICES" begin
		V = [[0.,0.] [1,0] [2,0] [0,1] [1,1] [2,1] [0,2] [1,2] [2,2]];
		FV = [[1,2,4],[2,3,5],[3,5,6],[4,5,7],[5,7,8],[6,8,9]];
		model = (V,FV);
		pattern = repeat([1,2,-3],outer=4);		
		W,FW = EXTRUDESIMPLICES(model, pattern);
		@test typeof(EXTRUDESIMPLICES(model, pattern))==
			Tuple{Matrix{Float64},Cells}
		@test typeof(W)==Matrix{Float64}
		@test size(W)==(3,117)
		@test typeof(FW)==Cells
		@test length(FW)==144
	end

	@testset "SIMPLEXGRID" begin
		@test SIMPLEXGRID([1])==([0.0 1.0], Array{Int64,1}[[1, 2]])
		@test SIMPLEXGRID([1,1])==([0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0], 
			Array{Int64,1}[[1, 2, 3], [2, 3, 4]])
		@test typeof(SIMPLEXGRID([1,1,1])[2])==Array{Array{Int64,1},1}
		@test length(SIMPLEXGRID([1,1,1])[2])==6
		V,CV = SIMPLEXGRID([10,10,1])
		@test typeof(V)==Array{Float64,2}
		@test typeof(CV)==Array{Array{Int64,1},1}
		@test size(V)==(3,242)
		@test length(CV)==600
	end

end



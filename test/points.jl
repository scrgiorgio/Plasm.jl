using Test
using Plasm


@testset "points.jl" begin

   @testset "meaning" begin
      A = collect(1.0:18.0)
      V = reshape(A, 2,:)
      @test size(V) == (2,9)
      V = reshape(A, 3,:)
      @test size(V) == (3,6)
      @test PDIM(V) == size(V,1)
      @test NVERTS(V) == size(V,2)
      @test size(BYROW(V)) ==  size(V')
      @test size(permutedims(V)) ==  size(V')
      @test BYROW(V) == permutedims(V)
   end;

   @testset "space" begin
      A = rand(Float64, 24)
      length(A) == 3*8 == 2*12
      V = reshape(A, 2,:)
      @test size(V) == (2,12)
      V = reshape(A, 3,:)
      @test size(V) == (3,8)
      @test PDIM(V) == size(V,1)
      @test NVERTS(V) == size(V,2)
      @test size(BYROW(V)) ==  size(V')
      @test size(permutedims(V)) ==  size(V')
      @test BYROW(V) == permutedims(V)
   end;

   @testset "fuzzy vertex" begin
      Plasm.LAR_DEFAULT_ERR,old_value = 1.0e-9, Plasm.LAR_DEFAULT_ERR
      v1 = v2 = rand(3)
      v1 = v2 .- 0.00000000009
      @test vertex_fuzzy_equals(v1,v2)
      w1 = w2 = rand(3);
      w1 = w2 .- 0.00000000009
      @test vertex_fuzzy_equals(w1,w2)
      t1 = t2 = rand(Float32,3);
      t1 = t2 .- 0.00000000009
      @test vertex_fuzzy_equals(t1,t2)
      Plasm.LAR_DEFAULT_ERR=old_value
   end;

   @testset "vertex membership in vertices" begin
      vertices = rand(0.000000001:0.00000000001:0.000000002, 10)
      vertex = rand(0.0000000001:0.00000000001:0.0000000002, 1)
      @test is_visited_vertex(vertex, vertices)
   end;

end;
using Plasm, Test
using SparseArrays
import Base.size

@testset "sparse.jl" begin

   obj = CUBOIDGRID([1,1,1])
   V,FV,EV = obj.V, obj.C[:FV], obj.C[:EV]
   
   @testset "lar2cop" begin
      copEV = lar2cop(EV)
      copFV = lar2cop(FV)
      @test typeof(copFV) == ChainOp == SparseMatrixCSC{Int8, Int64}
      @test SparseArrays.size(copFV) == (6,8)
      @test count(! iszero, copFV)  == 24
      @test typeof(lar2cop(FV)) == typeof(lar2cop(EV))
      @test typeof(copEV) == ChainOp == SparseMatrixCSC{Int8, Int64}
      @test SparseArrays.size(copEV) == (12,8)
      @test count(! iszero, copEV) == 24
      copFE = SparseArrays.spmatmul(copFV,permutedims(copEV)) .÷ Int8(2)
      @test SUM(sum(copFE, dims = 1)) == 24
      @test SUM(sum(copFE, dims = 2)) == 24
      SparseArrays.size(copFE) = (6,12)
      count(! iszero, copFE) == 24
      I = [1, 4, 3, 5]; J = [4, 7, 18, 9]; V = [1, 2, -5, 3];
      S = sparse(I,J,V)
      triplevec = findnz(S)
      @test findnz(S) == ([1, 4, 5, 3], [4, 7, 9, 18], [1, 2, 3, -5])
      @test findall(!iszero, S) == CartesianIndex{2}[CartesianIndex(1, 4),   CartesianIndex(4, 7), CartesianIndex(5, 9), CartesianIndex(3, 18)]
      @test triplevec[1]==[1,4,5,3]
      @test triplevec[2]==[4, 7, 9, 18]
      @test triplevec[3]==[1, 2, 3, -5]
   end;
   
   @testset "cop2lar" begin
      grid = CUBOIDGRID([2,2,1])
      V,FV,EV = grid.V, grid.C[:FV], grid.C[:EV]
      copFV, copEV = lar2cop(FV), lar2cop(EV)
      @test typeof(findnz(copFV)) == Tuple{Vector{Int64},Vector{Int64}, 
      Vector{Int8}}
      @test typeof(findnz(copEV)) == Tuple{Vector{Int64}, Vector{Int64}, 
      Vector{Int8}}
      @test cop2lar(copFV) == [[1, 2, 3, 4], [3, 4, 5, 6], [1, 3, 5, 7], [2, 4, 6, 8], [1, 2, 7, 8], [5, 6, 7, 8], [1, 3, 9, 10], [3, 5, 10, 11], [9, 10, 11, 12], [1, 7, 9, 12], [5, 7, 11, 12], [1, 2, 13, 14], [1, 7, 13, 15], [2, 8, 14, 16], [13, 14, 15, 16], [7, 8, 15, 16], [1, 9, 13, 17], [9, 12, 17, 18], [13, 15, 17, 18], [7, 12, 15, 18]]
      @test length(cop2lar(copFV)) == 20
      @test length(cop2lar(copEV)) == 33
      copFE = SparseArrays.spmatmul(copFV,permutedims(copEV)) .÷ Int8(2)
      @test length(cop2lar(copFE)) == 20
      FE = cop2lar(copFE)
      @test length(FE) == length(FV) == 20
   end;
   
   @testset "FV2EV" begin
      grid = CUBOIDGRID([2,2,1])
      V,FV,EV = grid.V, grid.C[:FV], grid.C[:EV]
      copFV, copEV = lar2cop(FV), lar2cop(EV)
      copFE = SparseArrays.spmatmul(copFV,permutedims(copEV)) .÷ Int8(2)
      ev = FV2EV(copEV::ChainOp, copFE::ChainOp)
      @test sort(ev) == sort(EV)
      @test typeof(ev) == Cells
      @test LEN(EV) == LEN(ev)
   end;
   
   @testset "find_vcycle" begin
      grid = CUBOIDGRID([2,2,1])
      V,FV,EV = grid.V, grid.C[:FV], grid.C[:EV]
      copFV, copEV = lar2cop(FV), lar2cop(EV)
      copFE = SparseArrays.spmatmul(copFV, permutedims(copEV)) .÷ Int8(2)
      copFV1 = SparseArrays.spmatmul(copFE, copEV) .÷ Int8(2)
      copFV2 = lar2cop(FV)
      @test copFV1 == copFV2 # calcolato con due diversi metodi
      vs, edges = find_vcycle(copEV::ChainOp, copFE::ChainOp, 1)
      @test LEN(vs) == LEN(edges)
      vs, edges = find_vcycle(copEV::ChainOp, copFE::ChainOp, copFE.m)
      @test LEN(vs) == LEN(edges)
   end;
   
end;
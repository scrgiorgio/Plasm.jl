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
   end;
   
   @testset "FV2EV" begin
   end;
   
   @testset "find_vcycle" begin
   end;
   
end;
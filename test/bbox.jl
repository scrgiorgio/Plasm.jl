using Plasm, Test
using IntervalTrees


#@testset "bbox.jl" begin

   obj = CUBOIDGRID([2,3])
   V,FV,EV = obj.V, obj.C[:FV], obj.C[:EV]
   
   @testset "lar_bbox_contains" begin
      @test size(V) == (2,12)
      @test PDIM(V) == size(V,1)
      @test NVERTS(V) == size(V,2)
      @test size(BYROW(V)) ==  size(V')
      @test size(permutedims(V)) ==  size(V')
      @test BYROW(V) == permutedims(V)
   end;

   @testset "lar_bbox_create" begin
      pmin, pmax = lar_bbox_create(V::Points; dims=2)
      @test typeof(pmin) == Matrix{Float64}
      @test typeof(pmax) == Matrix{Float64}
      W = V[:,FV[end]]
      fpmin,fpmax = lar_bbox_create(W::Points; dims=2)
      @test typeof(fpmin) == typeof(fpmax)
      @test size(fpmin) == size(fpmax)
      @test  all(fpmin .<= fpmax) == true
      @test !all(fpmin .>= fpmax) == true
   end;

   @testset "lar_bbox_contains-2" begin
      W = V[:,FV[1]]
      pmin, pmax = lar_bbox_create(V::Points; dims=2) # su tutto V
      fpmin,fpmax = lar_bbox_create(W::Points; dims=2) # sul sottoinsieme W
      b1_min, b1_max = pmin, pmax
      b2_min, b2_max = fpmin,fpmax
      @test all(map((i, j, k, l) -> i <= j <= k <= l, b1_min, b2_min, b2_max, 
      b1_max)) # a suo tempo ho fatto moltissimi test sulla coincidenza 
      W = V[:,EV[5]]
      emin, emax = lar_bbox_create(V::Points; dims=2)
      b2_min, b2_max = emin, emax
      @test all(map((i, j, k, l) -> i <= j <= k <= l, b1_min, b2_min, b2_max, 
      b1_max))
   end;

   @testset "lar_bbox_covering" begin
   # lo farò più tardi ...
   end;

   @testset "lar_bbox_coord_intervals" begin
   # lo farò più tardi ...
   end;

   @testset "lar_bbox_containment_graph" begin
   # ne vedo l'utilità solo nella grafica ma non qui
   # le ns operazioni booleane non sono basate su chi interseca chi ...
   # questo aspetto è già stato affrontato nello splitting 2D
   # non ha alcun interesse adesso ... !!!
   end;

#end;




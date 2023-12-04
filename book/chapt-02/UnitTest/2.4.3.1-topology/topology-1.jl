using LinearAlgebra

tetra::Lar = simplex(3) #=
Lar(3, 3, 4, [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0], Dict{Symbol, AbstractArray}(:c3v => [[1, 2, 3, 4]])) =#
typeof(tetra) # Lar
tetra.C #=
Dict{Symbol, AbstractArray} with 1 entry:
  :c3v => [[1, 2, 3, 4]] =#
tetra.V #=
3×4 Matrix{Float64}:
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0 =#
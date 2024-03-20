using Plasm
include("./arrangement.jl")
using SparseArrays

function CSC( CV ) # CV => Cells defined by their Vertices
   I = vcat( [ [k for h in CV[k]] for k=1:length(CV) ]...)                
   # vcat maps arrayofarrays to single array
   J = vcat( CV...)         
   # splat transforms the CV array elements to vcat arguments
   X = Int8[1 for k=1:length(I)]         
   # Type Int8 defines the memory map of array elements
   return SparseArrays.sparse(I,J,X)        
end

mutable struct Lar
  d::Int # intrinsic dimension
  m::Int # embedding dimension (rows of V)
  n::Int # number of vertices  (columns of V)
  V::Matrix{Float64} # object geometry
  C::Dict{Symbol, AbstractArray} # object topology (C for cells)
  # inner constructors
  Lar() = new( -1, 0, 0, Matrix{Float64}(undef,0,0), Dict{Symbol, AbstractArray}() )
  Lar(m::Int,n::Int) = new( m,m,n, Matrix(undef,m,n), Dict{Symbol,AbstractArray}() )
  Lar(d::Int,m::Int,n::Int) = new( d,m,n, Matrix(undef,m,n), Dict{Symbol,AbstractArray}() ) 
  Lar(V::Matrix) = begin m, n = size(V); 
     new( m,m,n, V, Dict{Symbol,AbstractArray}() ) end
  Lar(V::Matrix,C::Dict) = begin m,n = size(V); new( m,m,n, V, C )  end
  Lar(d::Int,V::Matrix,C::Dict) = begin m,n = size(V); new( d,m,n, V, C )  end
  Lar(d,m,n, V,C) = new( d,m,n, V,C )
end

function update(obj::Lar, (a,b))::Lar
   obj.C[a] = b
   return obj
end


esa = CUBE(1)
Lar(esa::Hpc)::Lar
twocubes = STRUCT([esa,T(3)(1),esa])
doublecube = Lar(twocubes::Hpc)::Lar
KEV = lar2cop(doublecube.C[:EV])

update(doublecube::Lar, (:KEV, KEV))
doublecube.C
doublecube.V

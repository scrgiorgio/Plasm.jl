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

""" Remove possible duplicates from ToLAR faces """
function removedups(obj::Hpc)::Cells
   # initializations
   geo=ToGeometry(obj)
   hulls = geo.faces
   dict  = geo.db
   inverted_dict = Dict{valtype(dict), Vector{keytype(dict)}}()
   [push!(get!(() -> valtype(inverted_dict)[], inverted_dict, v), k) for (k, v) in dict]  
   DB = []  # convert arrays of indices to arrays of points
   for hull in hulls
      points = []
      [ append!(points, inverted_dict[k]) for k in hull ]
      push!(DB, Set(points))
   end 
   DB = Set(DB) # eliminate duplicates
   faces = [[dict[point] for point in set] for set in DB]
   faces = sort!(AA(sort!)(faces))
end

""" Alternate implementation: Remove possible duplicates """
function removedups(obj::Hpc)::Cells
   # initializations
   geo=ToGeometry(obj)
   hulls = geo.faces
   dict  = geo.db
   inverted_dict = Dict{valtype(dict), Vector{keytype(dict)}}()
   [push!(get!(() -> valtype(inverted_dict)[], inverted_dict, v), k) for (k, v) ∈ dict]  
   # convert arrays of indices to arrays of points
   DB = [Set([inverted_dict[k] for k ∈ hull]) for hull ∈ hulls]
   DB = Set(DB) # eliminate duplicates
   faces = [[dict[point...] for point in set] for set ∈ DB]
   faces = sort!(AA(sort!)(faces))
end

function Lar(obj::Hpc)::Lar
   geo=ToGeometry(obj)
    V  = geo.points
    CV = geo.hulls
    faces = geo.faces
    FV = removedups(obj)
    FF = CSC(FV) * CSC(FV)'
    edges = filter(x->x[1]<x[2] && FF[x...]==2,collect(zip(findnz(FF)[1:2]...)))
    EV = sort!(collect(Set([FV[i] ∩ FV[j] for (i,j) in edges])))
    out = Lar(hcat(V...), Dict(:CV=>CV, :FV=>FV, :EV=>EV))
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

using Plasm
include("./arrangement.jl")

""" Remove possible duplicates from ToLAR facets """
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

""" Remove possible duplicates from ToLAR facets """
function removedups(obj::Hpc)::Cells
   # initializations
   hulls = ToLAR(obj).childs[1].facets
   dict = ToLAR(obj).childs[1].db
   inverted_dict = Dict{valtype(dict), Vector{keytype(dict)}}()
   for (k, v) in dict
      push!(get!(() -> valtype(inverted_dict)[], inverted_dict, v), k)
   end
   # convert arrays of indices to arrays of points
   DB = []
   for hull in hulls
      points = []
      for k in hull
          append!(points, inverted_dict[k])
      end
      push!(DB, Set(points))
   end
   # eliminate duplicates
   DB = Set(DB)
   # inverse conversion from arrays of points to arrays of indices
   faces = []
   for set in DB
      face = Int[]
      for point in set
         push!(face, dict[point])
      end
      push!(faces, face)
   end   
   faces = sort!(AA(sort!)(faces))
end

function Lar(obj::Hpc)::Lar
    V = ToLAR(obj).childs[1].points
    CV = ToLAR(obj).childs[1].hulls
    facets = ToLAR(obj).childs[1].facets
    FV = removedups(obj)
    FF = lar2cop(FV) * lar2cop(FV)'
    edges = filter(x->x[1]<x[2] && FF[x...]==2,collect(zip(findnz(FF)[1:2]...)))
    EV = sort!(collect(Set([intersect(FV[i],FV[j]) for (i,j) in edges])))
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

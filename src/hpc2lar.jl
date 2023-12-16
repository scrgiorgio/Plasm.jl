# //////////////////////////////////////////////////////////////////////////////
# ///Hpc -> Lar -> HPC  User Data type conversions//////////////////////////////
# //////////////////////////////////////////////////////////////////////////////
using DataStructures, SparseArrays

export truncate, simplifyCells, CSC, Lar, Hpc

# //////////////////////////////////////////////////////////////////////////////
# With a docstring, can be seen from Help (?) the whole set of parameters
"""
	truncate(PRECISION::Int)(value::Float64)

Transform the float `value` to get a `PRECISION` number of significant digits.
"""
truncate = PRECISION -> value -> begin
   approx = round(value,digits=PRECISION)
   abs(approx)==0.0 ? 0.0 : approx
end

# //////////////////////////////////////////////////////////////////////////////
"""
	W,CW = simplifyCells(V,CV)

Find and remove the duplicated vertices and the incorrect cells.

Some vertices could appear two or more times, due to numerical errors
on mapped coordinates. So, close vertices are identified, according to the
PRECISION number of significant digits.
"""
function simplifyCells(V,CV)
	PRECISION = 14
	vertDict = DefaultDict{Vector{Float64}, Int64}(0)
	index = 0
	W = Vector{Float64}[]
	FW = Vector{Int64}[]

	for incell in CV
		outcell = Int64[]
		for v in incell
			vert = V[:,v]
			key = map(truncate(PRECISION), vert)
			if vertDict[key]==0
				index += 1
				vertDict[key] = index
				push!(outcell, index)
				push!(W,key)
			else
				push!(outcell, vertDict[key])
			end
		end
		append!(FW, [[Set(outcell)...]])
	end
	return hcat(W...), filter(x->!(LEN(x)<3), FW)
end

# //////////////////////////////////////////////////////////////////////////////
"""
   CSC( Cc::Vector{Vector{Int64}} )::SparseMatrix

Creation of Compressed Sparse Column (CSC) sparse matrix format.
Each CV element is the array of vertex indices of a cell.
"""
function CSC( CV ) # CV => Cells defined by their Vertices
   I = vcat( [ [k for h in CV[k]] for k=1:length(CV) ]...)                
   # vcat maps arrayofarrays to single array
   J = vcat( CV...)         
   # splat transforms the CV array elements to vcat arguments
   X = Int8[1 for k=1:length(I)]         
   # Type Int8 defines the memory map of array elements
   return SparseArrays.sparse(I,J,X)        
end


# //////////////////////////////////////////////////////////////////////////////
"""
   mutable struct Lar

Linear Algebraic Representation (LAR). Data type for Cellular and Chain Complex.
"""
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

# //////////////////////////////////////////////////////////////////////////////
"""
   Lar(obj::Hpc)::Lar

Constructor of object of Linear Algebraic Representation (Lar) type, from Hpc object. 
"""
function Lar(obj::Hpc)::Lar
   V = ToLAR(obj).childs[1].points
   #CV = obj.childs[1].hulls
   facets = ToLAR(obj).childs[1].facets
   V,FV = simplifyCells( hcat(V...), facets )
   FF = CSC(FV) * CSC(FV)'
   edges = filter(x->x[1]<x[2] && FF[x...]==2,collect(zip(findnz(FF)[1:2]...)))
   EV = sort!(collect(Set([ FV[i] ∩ FV[j] for (i,j) in edges ]))) # ∩ intersect
   out = Lar( V, Dict(# :CV=>CV, 
      :FV=>FV, :EV=>EV) )
end


# //////////////////////////////////////////////////////////////////////////////
"""
   Hpc(V,CV)::Hpc 

Constructor of object of Hierarchical Polyhedral Complex (Hpc) type, starting from a pair V,CV of LAR kind. 
V is of type Matrix{Float64}; CV is any ::Vector{Vector{Int}} dataset.
"""
function Hpc(V::Matrix{Float64}, CV::Vector{Vector{Int}}) 
   W = [V[:,k] for k=1:size(V,2)]
   out = STRUCT(AA(MKPOL)(DISTL(W, AA(LIST)(CV)))) 
   return out 
end

function Hpc(obj::Lar) 
   V = obj.V;  FV = obj.C[:FV]
   return Hpc(V,FV) 
end


using LinearAlgebra

using DataStructures

mutable struct HPC
          d::Int # intrinsic dimension
          m::Int # embedding dimension (rows of V)
          n::Int # number of vertices  (columns of V)
          V::Matrix{Float64} # object geometry
          C::Dict{Symbol, AbstractArray} # object topology (C for cells)
          # inner constructors
          HPC() = new( NaN, NaN, NaN, Matrix{Float64}, Dict{Symbol, AbstractArray} )
          HPC(m,n) = new( m,m,n, Matrix(undef,m,n), Dict{Symbol,AbstractArray} )
          HPC(d,m,n) = begin V = Matrix(undef,m,n), C = Dict{Symbol,AbstractArray}(); 
                    new( d,m,n, V, C ) end
          HPC(V) = begin C = Dict{Symbol,AbstractArray}(); m, n = size(V); 
             new( m,m,n, V, C ) end
          HPC(V,C) = begin m,n = size(V); new( m,m,n, V, C )  end
          HPC(d,V,C) = begin m,n = size(V); new( d,m,n, V, C )  end
          HPC(d,m,n, V,C) = new( d,m,n, V,C )
          HPC() = new( -1,0,0, Matrix(undef,0,0),Dict() )
       end

function simplex(n; complex=false)
        V = [zeros(n,1) I]
        CV = [collect(1:n+1)]
        C = Dict(Symbol("c$(n)v") => CV)
        if complex == false return HPC(V,C)
        else
           cells = CV
           for k=n:-1:1
             cells = simplexfacets(cells)
                key = Symbol("c$(k-1)v")
                push!(C, (key => cells) )
           end
           return HPC(V,C)
        end
       end

function hpcprod( hpc1, hpc2 )
          key1 = findmax(keys, hpc1.C)[2]
          key2 = findmax(keys, hpc2.C)[2]
          (V, cells1) = hpc1.V, values(hpc1.C[key1])
          (W, cells2) = hpc2.V, values(hpc2.C[key2])
          vertices = vertprod(V, W)
          thecells = cellprod(cells1, cells2, V, W)
          cells = [[vertices[v] for v in cell] for cell in thecells]
          verts = hcat(keys(vertices)...); d = size(verts, 1)
          return HPC(verts, Dict(Symbol("c$(d)v") => cells))
       end

function vertprod(V, W)
          vertices = DataStructures.OrderedDict(); k = 0
          for j in 1:size(V,2)
             v = V[:,j]
             for i in 1:size(W,2)
                w = W[:,i]
                id = [v; w]; k += 1
                vertices[id] = k
             end
          end
          return vertices
       end

function cellprod(cells1, cells2, V, W)
          cells = Vector[]
          for c1 in cells1
             for c2 in cells2
                cell = Vector[]
                for vc in c1
                   for wc in c2
                      push!(cell, [V[:,vc];W[:,wc]] )
                   end
                end
                append!(cells, [cell])
             end
          end
          return cells
       end

line = simplex(1)

square = hpcprod(line, line)

cube = hpcprod(square, line)

doubleline = HPC([0 1 2], Dict(:c1v => [[1,2],[2,3]]))

doublesquare = hpcprod(doubleline,line)

doublecube = hpcprod(doublesquare,line)

import Base: *



using LinearAlgebraicRepresentation

Lar = LinearAlgebraicRepresentation
LinearAlgebraicRepresentation


function cubegrid( shape; skeltns=false )
   d = length(shape)
   vertgrid = imageverts(shape)
   gridmap = gridskeleton(shape)
   cells = !skeltns ? gridmap(d) : [gridmap(id) for id=1:d+1]
   return hcat(map(collect,vertgrid)...), cells
end

apply(fun,a) = fun(a)
columnlist(arr) = [arr[:,k]  for k in 1:size(arr,2)]

function gridskeleton(shape)
  n = length(shape)
  function gridskeleton1( d::Int )
    @assert d ≤ n
    components = filterbyorder(n)[d.+1]
    componentcellists = [[map(f,x) for (f,x) in zip([grid(dim-1)
      for dim in shape], convert(Array{Int64,1},component))]
        for component in components]
    out = ([cellprod(map(columnlist,cellists)) for cellists in componentcellists])
    return vcat(out...)
  end
  return gridskeleton1
end


function cellprod(cellists)
   shapes = [length(item) for item in cellists]
   subscripts = cart([collect(range(0, length=shape)) for shape in shapes])
   indices = [collect(tuple) .+ 1 for tuple in subscripts]

   jointCells = [cart([cells[k] for (k,cells) in zip(index,cellists)])
                                 for index in indices]
   convertIt = [ (length(cellists[k][1]) > 1) ? shape .+ 1 : shape
    for (k,shape) in enumerate(shapes) ]
   [vcat(map(convertIt, map(collect,jointCells[j]))...) for j in 1:length(jointCells)]
end

function imageverts( shape::Vector{Int64} )::Matrix{Int64}
   vertexdomain(n) = hcat([k for k in 0:n-1]...)
   vertlists = [vertexdomain(k+1) for k in shape]
   vertgrid = vertprod(vertlists)
   return vertGrid
end

function index2addr( shape::Array{Int64,2} )
  n = length(shape)
  theShape = append!(shape[2:end],1)
  weights = [prod(theShape[k:end]) for k in range(1, length=n)]
  function index2addr0( multiIndex::Array{Int,1} )::Int
      return LinearAlgebra.dot(collect(multiIndex), weights) + 1
  end
  return index2addr0
end


function filterbyorder(n::Int)Array{Array{Array{Int8,1},1},1}
    terms = [[parse(Int8,bit) for bit in collect(term)] for term in binaryRange(n)]
    return [[term for term in terms if sum(term) == k] for k in 0:n]
end

function grid(n::Int)
   d -> d==0 ? [i for i=0:n+1] : hcat([[i,i+1] for i=0:n]...)
end

function binaryRange(n)
    return string.(range(0, length=2^n), base=2, pad=n)
end

function grid(n::Int)
   d -> d==0 ? [i for i=0:n+1]' : hcat([[i,i+1] for i=0:n]...)
end

apply(fun,a) = fun(a)

columnlist(arr) = [arr[:,k]  for k in 1:size(arr,2)]

function cart(args)::Array{Tuple,1}
   return sort(vcat(collect(Iterators.product(args...))...))
end

function collist(arr)
   [arr[:,k]  for k in 1:size(arr,2)] 
end
function larGridSkeleton(shape)
    n = length(shape)
    function larGridSkeleton0( d::Int )
    	  @assert d<=n
        components = Lar.filterByOrder(n)[d .+ 1]
		  componentCellLists = [ [map(f,x)  for (f,x) in  zip( [Lar.larGrid(dim)
			 for dim in shape], convert(Array{Int64,1},component) ) ]
				for component in components ]
        out = [ Lar.larCellProd(map(colList,cellLists)) for cellLists in componentCellLists ]
        return vcat(out...)
    end
    return larGridSkeleton0
end

function gridskeleton(shape)
    n = length(shape)
    function gridskeleton0( d::Int )
        components = filterbyorder(n)[d.+1]
		  componentcelllists = [[map(f,x) for (f,x) in zip([grid(dim)
			 for dim in shape], convert(Array{Int64,1},component))]
				for component in components]
        out = [cellprod(map(collist,celllists)) for celllists in componentcelllists]
        return vcat(out...)
    end
    return gridskeleton0
end

for celllists in componentcelllists
   cellprod(map(collist,celllists))
end

celllists = componentcelllists


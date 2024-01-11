using LinearAlgebra
using QHull

export IsPolytope, IsSimplex, simplex, simplexfacets, CHULL, CUBOIDGRID, GRID1

import Base.+  
+(f::Function, g::Function) = (x...) -> f(x...) + g(x...)  

import Base.*
*(pol1::Hpc, pol2::Hpc) = Power(pol1, pol2)

# //////////////////////////////////////////////////////////////////////////////
IsPolytope = AND ∘ CONS([  # pure FL style
   ISPOL, 
   EQ ∘ CONS([LEN ∘ S2 ∘ UKPOL, K(1)])
])

# //////////////////////////////////////////////////////////////////////////////
IsSimplex = AND ∘ CONS([  # pure FL style
   IsPolytope, 
   EQ ∘ CONS([LEN ∘ S1 ∘ UKPOL, RN + K(1)]) 
])

# /////////////////////////////////////////////////////////////////////////////
function simplex(d; complex=false)
   V = [zeros(d,1) I]
   CV = [collect(1:d+1)]
   C = Dict(Symbol("C$(d)V") => CV)
   if complex == false return Lar(V,C)
   else
      cells = CV
      for k=d:-1:1
      	 cells = simplexfacets(cells)
	 	 key = Symbol("C$(k-1)V")
         push!(C, (key => cells) )
      end
      return Lar(V,C)
   end
end  

# //////////////////////////////////////////////////////////////////////////////
function simplexfacets(simplices)
   @assert hcat(simplices...) isa Matrix
   out = Array{Int64,1}[]
	for simplex in simplices
		for v in simplex
			facet = setdiff(simplex,v)
			push!(out, facet)
		end
	end
	# remove duplicate facets
	return sort(collect(Set(out)))
end

# //////////////////////////////////////////////////////////////////////////////
function CHULL(points::Matrix)
   convexhull = QHull.chull(points)
   FV = convexhull.simplices
   pairs = [map(sort,[[i,j],[i,k],[j,k]]) for (i,j,k) in FV]
   EV = sort!(union(vcat(pairs...)))
   ret = Lar(2, Matrix(p'), Dict(:EV => EV, :FV => FV))
end

# //////////////////////////////////////////////////////////////////////////////
GRID1(n) = QUOTE(DIESIS(n)(1.0))

# //////////////////////////////////////////////////////////////////////////////
function CUBOIDGRID(shape::Vector{Int})
   LAR(INSL(POWER)(AA(GRID1)(shape)))
end

const Cells = Vector{Vector{Int}}
export Cells

# Linear Algebraic Representation . Data type for Cellular and Chain Complex.
mutable struct Lar

	# object geometry
	# always stored by column i.e. a point is a column 
	#     x1 x2 x3 ...
	#     y1 y2 y3 ...
	V::Points

	# object topology (C for cells)
	C::Dict{Symbol,AbstractArray}

	# constructor
	Lar(V::Matrix{Float64}=Matrix{Float64}(undef, 0, 0), C::Dict=Dict{Symbol,AbstractArray}()) = begin
		new(V, C)
	end

end
export Lar


# //////////////////////////////////////////////////////////////////////////////
"""from Hpc -> Lar 
"""
function LAR(obj::Hpc; precision=DEFAULT_PRECISION)::Lar
	geo = ToGeometry(obj, precision=precision)
	n = length(geo.points)    # number of vertices  (columns of V)
	m = length(geo.points[1]) # embedding dimension (rows of V) i.e. number of coordinates
	ret = Lar()
	ret.V = hcat(geo.points...)
	ret.C[:EV] = geo.edges
	ret.C[:FV] = geo.faces
	ret.C[:CV] = geo.hulls
	return ret
end
export LAR

# //////////////////////////////////////////////////////////////////////////////
""" Create an Hpc from Lar 

use MKPOLS to specify what Hpc you want to build (like only edges or only 2d-faces)
"""
function HPC(lar::Lar)::Hpc

	if :FV in keys(lar.C) && length(lar.C[:FV])
		return MKPOLS(lar.V, lar.C[:FV])

	elseif :EV in keys(lar.C) && length(lar.C[:EV])
		return MKPOLS(lar.V, lar.C[:EV])

	else
		error("Empty Lar")
	end

end
export HPC

# //////////////////////////////////////////////////////////////////////////////
function CUBOIDGRID(shape::Vector{Int})::Lar
	obj = INSL(POWER)(AA(GRID1)(shape))
	if RN(obj) == 2
		geo = ToGeometry(obj)
		V,FV,EV = geo.points,geo.hulls,geo.edges
		return Lar(hcat(V...), Dict(:FV => FV, :EV => EV))
	else
		return LAR(obj)
	end
end
export CUBOIDGRID


# ////////////////////////////////////////////////////////////////
""" create LAR_SIMPLEX

see also fenv SIMPLEX function
"""
function LAR_SIMPLEX(d; complex=false)

	function simplexfacets(simplices)
		@assert hcat(simplices...) isa Matrix
		out = Array{Int64,1}[]
		for it in simplices
			for v in it
				facet = setdiff(it, v)
				push!(out, facet)
			end
		end
		# remove duplicate faces
		return sort(union(out))
	end

	V = [zeros(d, 1) I]
	CV = [collect(1:d+1)]
	C = Dict(Symbol("C$(d)V") => CV)
	if complex == false
		return Lar(V, C)
	else
		cells = CV
		for k = d:-1:1
			cells = simplexfacets(cells)
			key = Symbol("C$(k-1)V")
			push!(C, (key => cells))
		end
		return Lar(V, C)
	end
end
export LAR_SIMPLEX


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


# /////////////////////////////////////////////////////////////////////
function EXPLODECELLS(V::Points, cells::Vector{Cells}; scale_factor::Float64=1.2)::Vector{Hpc}
	
	ret = []
	for cell in cells

		cell_indices = sort(union(cell...))

		cell_vertices=V[:, cell_indices]

		p1=compute_centroid(cell_vertices)
		p2 = PDIM(V) == 2 ? 
			p1 .* [scale_factor; scale_factor] : 
			p1 .* [scale_factor; scale_factor; scale_factor]

		cell_vertices = cell_vertices .+ (p2 - p1)

		vdict = Dict(zip(cell_indices, 1:length(cell_indices)))
		hulls=[[vdict[v_index] for v_index in face] for face in cell]
		push!(ret, MKPOL(cell_vertices,hulls))
	end
	return ret
end
export EXPLODECELLS


# /////////////////////////////////////////////////////////
# scrgiorgio: NOT (!!!) well tested
function compute_VE(EV::Cells)

	num_vertices=maximum([maximum([v1,v2]) for (v1,v2) in EV])
  ret=convert(Cells,[[] for k in 1:num_vertices])
  for (E,(v1,v2)) in enumerate(EV)
    push!(ret[v1],E)
    push!(ret[v2],E)
  end
  return ret
end
export compute_VE

# /////////////////////////////////////////////////////////
# scrgiorgio: NOT (!!!) well tested
function compute_FE(EV::Cells, FV::Cells; double_check=false)
  VE=compute_VE(EV)
  ret = convert(Cells,[])
  for face in FV
    fe=[]
    for vertex_index in face
      for edge_index in VE[vertex_index]
        a,b=EV[edge_index]
        if a in face && b in face
          push!(fe,edge_index)
        end
      end
    end
    push!(ret,sort(union(fe)))
  end

  if double_check
    for (F,face) in enumerate(ret)
      for edge in face
        a,b=EV[edge]
        @assert(a in FV[F])
        @assert(b in FV[F])
      end
    end
  end

  return ret
end
export compute_FE

# //////////////////////////////////////////////////////////////////////////////
""" return ordered vertices  and edges of the 1-cycle f """
function find_vcycle_v2(EV::Cells, FE::Cells, f::Int)
	vpairs = [EV[e] for e in FE[f]]
	ordered = []
	(A, B), todo = vpairs[1], vpairs[2:end]
	push!(ordered, A)
	while length(todo) > 0
		found = false
		for (I, (a, b)) in enumerate(todo)
			if a == B || b == B
				push!(ordered, B)
				B = (b == B) ? a : b
				found = true
				deleteat!(todo, I)
				break
			end
		end
		@assert found
	end
	push!(ordered, ordered[1])
	edges = [[a, b] for (a, b) in zip(ordered[1:end-1], ordered[2:end])]
	return Array{Int}(ordered[1:end-1]), edges
end
export find_vcycle_v2

# //////////////////////////////////////////////////////////////////////////////
function RandomLine(size_min::Float64,size_max::Float64)
	size = size_min+rand()*(size_max-size_min)
	return STRUCT(
		T(1,2)(rand(2)...), 
		S([1,2])([size,size]), 
		R([1,2])(2*pi*rand()),
		Plasm.SQUARE(1)
	)
end
export RandomLine

# //////////////////////////////////////////////////////////////////////////////
function RandomBubble()
  vs = rand()
  vt = rand(2)
  return STRUCT(
    T(1,2)(vt...),
    S([1,2])([0.25*vs, 0.25*vs]), 
    CIRCUMFERENCE(1)(rand(3:32))
  )
end
export RandomBubble

# //////////////////////////////////////////////////////////////////////////////
function RandomCube(size_min::Float64,size_max::Float64)
  size = size_min+rand()*(size_max-size_min)
  return STRUCT(
    T(1,2,3)(rand(3)...), 
    S([1,2,3])([size,size,size]), 
    R([1,2])(2*pi*rand()),
    R([2,3])(2*pi*rand()),
    R([1,3])(2*pi*rand()),
    Plasm.CUBE(1) 
  )
end
export RandomCube

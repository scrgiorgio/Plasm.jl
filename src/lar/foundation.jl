using LinearAlgebra
using SparseArrays
using DataStructures
using StaticArrays
using IntervalTrees
using Triangulate

const Points = Matrix{Float64}
export Points



# //////////////////////////////////////////////////////////////////////////////
""" Predicate to check equality of two vertices (only used above) """
function vertex_fuzzy_equals(v1, v2;err = 10e-8)
	return length(v1) == length(v2) && all(map((x1, x2) -> -err < x1 - x2 < err, v1, v2))
end
export vertex_fuzzy_equals

""" Predicate to check membership of `vertex` in `vertices_set` array"""
function is_visited_vertex(vertex, vertices_set)::Bool
	for v in vertices_set
		if vertex_fuzzy_equals(vertex, v)
			return true
		end
	end
	return false
end
export is_visited_vertex


const Cells = Vector{Vector{Int}}
export Cells


# //////////////////////////////////////////////////////////////////////////////
"""
bbox(vertices::Points)

The axis aligned *bounding box* of the provided Matrix of n-dim `vertices`.
The box is returned as the pair of `Points` of two opposite corners.
"""
function bbox(vertices::Points)
	minimum = mapslices(x -> min(x...), vertices, dims=1)
	maximum = mapslices(x -> max(x...), vertices, dims=1)
	minimum, maximum
end
export bbox

function boundingbox(vertices::Points)
	minimum = mapslices(x -> min(x...), vertices, dims=2)
	maximum = mapslices(x -> max(x...), vertices, dims=2)
	return minimum, maximum
end
export boundingbox

function bbox_contains(container, contained)
	b1_min, b1_max = container
	b2_min, b2_max = contained
	all(map((i, j, k, l) -> i <= j <= k <= l, b1_min, b2_min, b2_max, b1_max))
end
export bbox_contains

function boxcovering(bboxes, index, tree)
	covers = [[] for k = 1:length(bboxes)]
	for (i, boundingbox) in enumerate(bboxes)
		extent = bboxes[i][index, :]
		iterator = IntervalTrees.intersect(tree, tuple(extent...))
		for x in iterator
			append!(covers[i], x.value)
		end
	end
	return covers
end
export boxcovering

""" Make dictionary of 1D boxes for IntervalTrees construction """
function coordintervals(coord, bboxes)
	boxdict = OrderedDict{Array{Float64,1},Array{Int64,1}}()
	for (h, box) in enumerate(bboxes)
		key = box[coord, :]
		if haskey(boxdict, key) == false
			boxdict[key] = [h]
		else
			push!(boxdict[key], h)
		end
	end
	return boxdict
end
export coordintervals



# ///////////////////////////////////////////////////////////////////////////
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


""" returns point dim """
function PDIM(lar::Lar)
	return size(lar.V,1)
end
export PDIM

""" returns number of Lar vertices """
function NVERS(lar::Lar)
	return size(lar.V,2)
end
export NVERS

# convert a vertex matrix from by-col (LAR default) to by-col
# assumption: input is by-col
function BYROW(V::Matrix{Float64})::Matrix{Float64}
	return permutedims(V)
end
export BYROW

# convert a vertex matrix from by-row to by-col (LAR default)
# assumption: input is by-row
function BYCOL(V::Matrix{Float64})::Matrix{Float64}
	return permutedims(V)
end
export BYCOL


# //////////////////////////////////////////////////////////////////////////////
# sparse representation

const Cell         = SparseVector{Int8,Int}
const Chain        = SparseVector{Int8,Int}
const ChainOp      = SparseMatrixCSC{Int8,Int}
const ChainComplex = Vector{ChainOp}
export Cell, Chain, ChainOp, ChainComplex


function characteristicMatrix(FV::Cells)::ChainOp
	I, J, V = Int64[], Int64[], Int8[]
	for f = 1:length(FV)
		for k in FV[f]
			push!(I, f)
			push!(J, k)
			push!(V, 1)
		end
	end
	return sparse(I, J, V)
end
export characteristicMatrix

""" converte dense represantation to sparse"""
function lar2cop(CV::Cells)::ChainOp
	I = Int[]
	J = Int[]
	Value = Int8[]
	for k = 1:size(CV, 1)
		n = length(CV[k])
		append!(I, k * ones(Int, n))
		append!(J, CV[k])
		append!(Value, ones(Int, n))
	end
	return SparseArrays.sparse(I, J, Value)
end
export lar2cop

""" converte sparse represantation to dense"""
function cop2lar(cop::ChainOp)::Cells
	[findnz(cop[k, :])[1] for k = 1:size(cop, 1)]
end
export cop2lar

function boundary_1(EV::Cells)::ChainOp
	out = characteristicMatrix(EV)'
	for e = 1:length(EV)
		out[EV[e][1], e] = -1
	end
	return out
end
export boundary_1

function coboundary_0(EV::Cells)::ChainOp
	return convert(ChainOp, LinearAlgebra.transpose(boundary_1(EV::Cells)))
end
export coboundary_0

# //////////////////////////////////////////////////////////////////////////////
"""Find EV from EV FE"""
function FV2EVs(copEV::ChainOp, copFE::ChainOp)
	EV = [findnz(copEV[k, :])[1] for k = 1:size(copEV, 1)]
	FE = [findnz(copFE[k, :])[1] for k = 1:size(copFE, 1)]
	EVs = [[EV[e] for e in fe] for fe in FE]
	return EVs
end
export FV2EVs


# //////////////////////////////////////////////////////////////////////////////
# from Hpc -> Lar (trying to keep as many information as possible, still unifying to a unique geometry)
export LAR
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

# //////////////////////////////////////////////////////////////////////////////
#NOTE: better use MKPOLS to specify what Hpc you want to build 
export HPC
function HPC(lar::Lar)::Hpc

	if :FV in keys(lar.C) && length(lar.C[:FV])
		return MKPOLS(lar.V, lar.C[:FV])

	elseif :EV in keys(lar.C) && length(lar.C[:EV])
		return MKPOLS(lar.V, lar.C[:EV])

	else
		error("Empty Lar")
	end

end

# //////////////////////////////////////////////////////////////////////////////
export LAR_CUBOIDGRID
function LAR_CUBOIDGRID(shape::Vector{Int})::Lar
	obj = INSL(POWER)(AA(GRID1)(shape))
	if RN(obj) == 2
		geo = ToGeometry(obj)
		V = geo.points
		FV = geo.hulls
		EV = geo.edges
		return Lar(hcat(V...), Dict(:FV => FV, :EV => EV))
	else
		return LAR(obj)
	end
end


# ////////////////////////////////////////////////////////////////
export LAR_SIMPLEX
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



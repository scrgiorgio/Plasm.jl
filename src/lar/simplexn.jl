#///////////////////////////////////////////////////////////////////////////////
"""
	SIMPLEX(n::Int; boundary=false::Bool)

Return a `LAR` model of the *`n`-dimensional simplex* in *`n`-space*.
When `fullmodel==true` return a `LARmodel`, including the faces, from dimension `1` to `n`.
# Example
```
V, cells = LARSIMPLEX(3,boundary=true)
```
"""
function LARSIMPLEX(n; boundary=false)
	eye(n) = LinearAlgebra.Matrix{Int}(I,n,n)
	V = [zeros(n,1) eye(n)]
	CV = [collect(1:n+1)]
	if boundary == false
		return V,CV
	else
		h = n
		cells = [CV]
		while h != 0
			push!(cells, SIMPLEXFACETS(cells[end]))
			h -= 1
		end
		return V,reverse(cells)
	end
end

export LARSIMPLEX

SIMPLEX=LARSIMPLEX
export SIMPLEX


#///////////////////////////////////////////////////////////////////////////////
"""
	EXTRUDESIMPLICES(model::Pair, pattern::Array)

Algorithm for multimensional extrusion of a simplicial complex.
Can be applied to 0-, 1-, 2-, ... simplicial models, to get a 1-, 2-, 3-, .... model.
The pattern::Array is used to specify how to decompose the added dimension.

A `model` is a LAR model, i.e. a pair (vertices,cells) to be extruded, whereas pattern is an array of `Int64`, to be used as lateral measures of the *extruded* model. `pattern` elements are assumed as either *solid* or *empty* measures, according to their (+/-) sign.

# Example
```julia
julia> V = [[0.,0] [1,0] [2,0] [0,1] [1,1] [2,1] [0,2] [1,2] [2,2]];
julia> FV = [[1,2,4],[2,3,5],[3,5,6],[4,5,7],[5,7,8],[6,8,9]];
julia> pattern = repeat([1,.2,-2],outer=4);
julia> model = (V,FV);
julia> W,FW = EXTRUDESIMPLICES(model, pattern);
julia> VIEW(MKPOL(W,FW))
julia> VIEWCOMPLEX(LAR(MKPOL(W,FW)))
```
"""
function EXTRUDESIMPLICES(model, pattern)
	 V = [model[1][:,k] for k=1:size(model[1],2)]
    FV = model[2]
    d, m = length(FV[1]), length(pattern)
    coords = collect(cumsum(append!([0.], abs.(pattern))))
    offset, outcells, rangelimit, i = length(V), [], d*m, 0
    for cell in FV
    	i += 1
        tube = [v+k*offset for k in range(0, length=m+1) for v in cell]
        cellTube = [tube[k:k+d] for k in range(1, length=rangelimit)]
        if i==1 outcells = reshape(cellTube, d, m)
        else outcells = vcat(outcells, reshape(cellTube, d, m)) end
    end
    cellGroups = []
    for i in 1:size(outcells, 2)
        if pattern[i]>0
            cellGroups = vcat(cellGroups, outcells[:, i])
        end
    end
    outVertices = [vcat(v, [z]) for z in coords for v in V]
    cellGroups = convert(Array{Array{Int, 1}, 1}, cellGroups)
    outModel = convert(PointsNd,outVertices), cellGroups
    hcat(outVertices...), cellGroups
end

export EXTRUDESIMPLICES

#///////////////////////////////////////////////////////////////////////////////
"""
	SIMPLEXGRID(shape::Array)::LAR

Generate a simplicial complex decomposition of a cubical grid of ``d``-cuboids, where ``d`` is the length of `shape` array. Vertices (0-cells) of the grid have `Int64` coordinates.

# Examples
```julia
julia> SIMPLEXGRID([0]) # 0-dimensional complex
julia> V,EV = SIMPLEXGRID([1]) # 1-dimensional complex
julia> V,FV = SIMPLEXGRID([1,1]) # 2-dimensional complex
julia> V,CV = SIMPLEXGRID([10,10,1]) # 3-dimensional complex
julia> VIEW(MKPOLS(V,CV))
julia> VIEWCOMPLEX(LAR(MKPOLS(V,CV)))
"""
function SIMPLEXGRID(shape)
    model = [], [[1]]
    for item in shape
        model = EXTRUDESIMPLICES(model, fill(1, item))
    end
    V, CV = model
    V = convert(Matrix{Float64}, V)
    CV = convert(Cells, CV)
    return V, CV
end

export SIMPLEXGRID


#///////////////////////////////////////////////////////////////////////////////

"""
	SIMPLEXFACETS(simplices::Cells)::Cells

Compute the `(d-1)`-skeleton (unoriented set of `facets`) of a simplicial `d`-complex.

# Example
```julia
julia> V,FV = SIMPLEXGRID([1,1]) # 2-dimensional complex
julia> W,CW = EXTRUDESIMPLICES((V,FV), [1,.5,1])
julia> FW = SIMPLEXFACETS(CW)
julia> V,CV = SIMPLEXGRID([3,3,3])
julia> VIEW(MKPOLS(V,CV))
julia> V,CV = SIMPLEXGRID([1,1,1])
julia> FV = SIMPLEXFACETS(CV)
julia> EV = SIMPLEXFACETS(FV)
julia> V = SIMPLEXFACETS(EV)
```
"""
function SIMPLEXFACETS(simplices)
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

export SIMPLEXFACETS

#///////////////////////////////////////////////////////////////////////////////
"""
	QUADS2TRIANGLES(quads::Cells)::Cells

Convert an array of *quads* with type `::Lar.Cells` into an array of *triangles*
with the same type.
# Examples
The transformation from quads to triangles works for any 2-complex, embedded in any dimensional space
```
julia> obj = Plasm.CUBOIDGRID([4,5])::Lar;
julia> triangles = QUADS2TRIANGLES(obj.C[:FV]::Cells)::Cells
julia> obj = Plasm.CUBOIDGRID([4,5,3])
julia> triangles = QUADS2TRIANGLES(obj.C[:FV]::Cells)::Cells
julia> VIEWCOMPLEX(LAR(MKPOLS(obj.V,triangles)))
julia> quads = obj.C[:FV]
julia> VIEWCOMPLEX(LAR(MKPOLS(obj.V,quads)))
```
"""
function QUADS2TRIANGLES(quads::Plasm.Cells)::Plasm.Cells
	pairs = [[ Int[v1,v2,v3], Int[v3,v2,v4]] for (v1,v2,v3,v4) in quads ]
	return CAT(pairs)
end

export QUADS2TRIANGLES
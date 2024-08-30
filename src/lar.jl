
# Linear Algebraic Representation . Data type for Cellular and Chain Complex.
export Lar
mutable struct Lar

	# object geometry
	# always stored by column i.e. a point is a column and
	#     x1 x2 x3 ...
	#     y1 y2 y3 ...
	V::Matrix{Float64}

	# object topology (C for cells)
	C::Dict{Symbol,AbstractArray}

	# constructor
	Lar(V::Matrix{Float64}=Matrix{Float64}(undef, 0, 0), C::Dict=Dict{Symbol,AbstractArray}()) = begin
		new(V, C)
	end

end

function PDIM(lar::Lar)
	return size(lar.V,1)
end

function NVERS(lar::Lar)
	return size(lar.V,2)
end

# convert a vertex matrix from by-col (LAR default) to by-col
export BYROW
function BYROW(V::Matrix{Float64})::Matrix{Float64}
	return permutedims(V)
end

# convert a vertex matrix from by-row to by-col (LAR default)
export BYCOL
function BYCOL(V::Matrix{Float64})::Matrix{Float64}
	return permutedims(V)
end

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
export LAR2TRIANGLES
""" input old LAR consistent data; output triangulated_faces """
function LAR2TRIANGLES(V, EV, FV, FE)

	""" return ordered vertices  and edges of the 1-cycle f """
	function __find_cycle(EV, FE, f::Int)
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

	triangulated_faces = Vector{Any}(undef, length(FE))

	for edges_idxs in FE
		edge_num = length(edges_idxs)
		fv, edges = __find_cycle(EV, FE, f)
		# look for independent vector triple
		points = V[:, fv]
		vmap = Dict(zip(fv, 1:length(fv))) # vertex map
		mapv = Dict(zip(1:length(fv), fv)) # inverse vertex map
		edges = [[vmap[A], vmap[B]] for (A, B) in edges]
		v1 = LinearAlgebra.normalize(points[2, :] - points[1, :])
		v2 = [0, 0, 0]
		v3 = [0, 0, 0]
		err = 1e-8
		i = 3
		while -err < LinearAlgebra.norm(v3) < err
			v2 = LinearAlgebra.normalize(points[i, :] - points[1, :])
			v3 = LinearAlgebra.cross(v1, v2)
			i = i % size(points, 1) + 1
		end
		# independent vector triple in face f 
		M = [v1 v2 v3]
		projected = BYCOL((points*M)[:, 1:2])
		trias = constrained_triangulation2D(projected, edges)  # single face f
		triangulated_faces[f] = [[mapv[v] for v in tria] for tria in trias]

	end
	return triangulated_faces
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

# //////////////////////////////////////////////////////////////////////////////
export POPOUTER
function POPOUTER(V, pols)

	# extract bounding boxes to find outser space
	bboxes = []
	for pol in pols
		ev = pol[1] # atom edges (pairs of vertex indices)
		verts = sort(union(CAT(CAT(ev)))) # atom vertices (vertex indices)
		v = V[:, verts]
		xbox = bbox(v) # maxmin of first coordinate
		push!(bboxes, xbox)
	end

	# compute 3D cuboidal volumes as (pmin,pmax)
	boxes = AA(collect ∘ AA(vec))(bboxes)
	diags = [LinearAlgebra.norm(v2 - v1) for (v1, v2) in boxes]
	__value, outer_position = findmax(diags)
	outerspace = filter(x -> x == pols[outer_position], pols)[1]
	atoms = filter(x -> x != pols[outer_position], pols)

	return outerspace, atoms
end

# ///////////////////////////////////////////////////
export VIEWCOMPLEX
function VIEWCOMPLEX(
	V,
	EV::Vector{Vector{Int64}},
	FV::Vector{Vector{Int64}};
	V_text::Vector{String}=nothing,
	EV_text::Vector{String}=nothing,
	FV_text::Vector{String}=nothing,
	properties::Properties=Properties())

	properties["background_color"] = get(properties, "background_color", DEFAULT_LAR_BACKGROUND_COLOR)
	properties["line_color"] = get(properties, "line_color", DEFAULT_LINE_COLOR)
	properties["line_width"] = get(properties, "line_width", DEFAULT_LINE_WIDTH)

	properties["text_v_color"] = get(properties, "text_v_color", DEFAULT_TEXT_V_COLOR)
	properties["text_ev_color"] = get(properties, "text_ev_color", DEFAULT_TEXT_EV_COLOR)
	properties["text_fv_color"] = get(properties, "text_fv_color", DEFAULT_TEXT_FV_COLOR)

	# font sizes
	properties["v_fontsize"] = get(properties, "v_fontsize", DEFAULT_V_FONTSIZE)
	properties["ev_fontsize"] = get(properties, "ev_fontsize", DEFAULT_EV_FONTSIZE)
	properties["fv_fontsize"] = get(properties, "fv_fontsize", DEFAULT_FV_FONTSIZE)

	if isnothing(V_text)
		V_text = [string(I) for I in eachindex(V)]
	end
	if isnothing(EV_text)
		EV_text = [string(I) for I in eachindex(EV)]
	end
	if isnothing(FV_text)
		FV_text = [string(I) for I in eachindex(FV)]
	end

	used_vertices = []
	for v in EV
		append!(used_vertices, v)
	end
	for v in FV
		append!(used_vertices, v)
	end

	obj = MKPOLS(V, EV)
	batches = Vector{GLBatch}()
	append!(batches, GetBatchesForHpc(obj))

	num_vertices = size(V, 2)
	vertices = [V[:, I] for I = 1:num_vertices]

	# show vertices
	if properties["text_v_color"][4] > 0.0 && properties["v_fontsize"] > 0
		for I in eachindex(vertices)
			append!(batches, GLText(
				(I in used_vertices ? V_text[I] : ""),
				center=ComputeCentroid([vertices[it] for it in [I]]),
				color=properties["text_v_color"],
				fontsize=properties["v_fontsize"]))
		end
	end

	# show edges
	if !isnothing(EV) && properties["text_ev_color"][4] > 0.0 && properties["ev_fontsize"] > 0
		for I in eachindex(EV)
			append!(batches, GLText(EV_text[I],
				center=ComputeCentroid([vertices[it] for it in EV[I]]),
				color=properties["text_ev_color"],
				fontsize=properties["ev_fontsize"]))
		end
	end

	# show faces
	if !isnothing(FV) && properties["text_fv_color"][4] > 0.0 && properties["fv_fontsize"] > 0
		for I in eachindex(FV)
			append!(batches,
				GLText(FV_text[I],
					center=ComputeCentroid([vertices[it] for it in FV[I]]),
					color=properties["text_fv_color"],
					fontsize=properties["fv_fontsize"]))
		end
	end

	View(batches, properties)

end

# //////////////////////////////////////////////////////////////////////////////
function VIEWCOMPLEX(mesh::Lar; properties::Properties=Properties())
	return VIEWCOMPLEX(
		mesh.V,
		:EV in keys(mesh.C) ? mesh.C[:EV] : nothing,
		:FV in keys(mesh.C) ? mesh.C[:FV] : nothing,
		properties=properties
	)
end

# //////////////////////////////////////////////////////////////////////////////
export SHOWEXPLODED
function SHOWEXPLODED(V, CVs, FVs, EVs; sx=1.2, sy=1.2, sz=1.2)

	# faces
	exploded = explodecells(V, FVs, sx=sx, sy=sx, sz=sx)
	v = []
	for k in 1:length(exploded)
		face_color = Point4d(Plasm.COLORS[(k-1)%12+1] - (rand(Float64, 4) * 0.1))
		face_color[4] = 1.0
		push!(v, PROPERTIES(exploded[k], Properties(
			"face_color" => face_color,
			#"line_color" => GL.BLACK, 
			"line_width" => 3)))
	end
	VIEW(STRUCT(v))

	# edges
	exploded = explodecells(V, EVs, sx=sx, sy=sx, sz=sx)
	v = []
	for k in 1:length(exploded)
		line_color = Point4d(Plasm.COLORS[(k-1)%12+1] - (rand(Float64, 4) * 0.1))
		line_color[4] = 1.0
		push!(v, PROPERTIES(exploded[k],
			Properties("line_color" => line_color, "line_width" => 3)))
	end
	VIEW(STRUCT(v))

	# full-dims 
	exploded = explodecells(V, CVs[1:end], sx=sx, sy=sx, sz=sx)
	v = []
	for k in 1:length(exploded)
		face_color = Point4d(Plasm.COLORS[(k-1)%12+1] - (rand(Float64, 4) * 0.1))
		face_color[4] = 1.0
		push!(v, PROPERTIES(exploded[k], Properties("line_width" => 3)))
	end
	VIEW(STRUCT(v))
end


# //////////////////////////////////////////////////////////////////////////////
export VIEWATOMS
function VIEWATOMS(V, copEV, copFE, copCF, atoms; view_outer=true)

	# per row
	num_vertices = size(V, 2)

	EV = AA(sort)(cop2lar(copEV))
	FE = AA(sort)(cop2lar(copFE))
	FV = AA(sort)([union(CAT([EV[e] for e in f])) for f in FE])

	ev_text = Dict{Vector{Int},Int}(collect(zip(EV, 1:length(EV))))
	fv_text = Dict{Vector{Int},Int}(collect(zip(FV, 1:length(FV))))

	outerspace, __others = POPOUTER(V, atoms)

	# Visualization of all atoms
	for atom in atoms

		# skip outer if needed
		if !view_outer && atom == outerspace
			continue
		end

		atom_EV, atom_FV = atom
		atom_EV = [sort(it) for it in union(CAT(atom_EV))]
		atom_FV = [sort(it) for it in atom_FV]

		VIEWCOMPLEX(
			V,
			atom_EV,
			atom_FV,
			V_text=[string(I) for I in 1:num_vertices],
			EV_text=[string(ev_text[it]) for it in atom_EV],
			FV_text=[string(fv_text[it]) for it in atom_FV]
		)
	end

end
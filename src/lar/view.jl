
# ///////////////////////////////////////////////////
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
export VIEWCOMPLEX

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
function VIEWEXPLODED(V, CVs, FVs, EVs; scale_factor=1.2)

	# faces
	v = []
	for (k,it) in enumerate(EXPLODECELLS(V, FVs, scale_factor=scale_factor))
		face_color = Point4d(Plasm.COLORS[(k-1)%12+1] - (rand(Float64, 4) * 0.1))
		face_color[4] = 1.0
		push!(v, PROPERTIES(it, Properties( "face_color" => face_color, "line_width" => 3)))
	end
	VIEW(STRUCT(v))

	# edges
	begin
		v = []
		for (k,it) in enumerate(EXPLODECELLS(V, EVs, scale_factor=scale_factor))
			line_color = Point4d(Plasm.COLORS[(k-1)%12+1] - (rand(Float64, 4) * 0.1))
			line_color[4] = 1.0
			push!(v, PROPERTIES(it, Properties("line_color" => line_color, "line_width" => 3)))
		end
		VIEW(STRUCT(v))
	end

	# full-dims 
	begin
		v = []
		for (k,it) in enumerate(EXPLODECELLS(V, CVs[1:end], scale_factor=scale_factor))
			face_color = Point4d(Plasm.COLORS[(k-1)%12+1] - (rand(Float64, 4) * 0.1))
			face_color[4] = 1.0
			push!(v, PROPERTIES(it, Properties("line_width" => 3, "face_color" => face_color)))
		end
		VIEW(STRUCT(v))
	end
end
export VIEWEXPLODED

# //////////////////////////////////////////////////////////////////////////////
function VIEWATOMS(V, copEV, copFE, copCF, atoms; view_outer=true)

	# per row
	num_vertices = size(V, 2)

	EV = AA(sort)(cop2lar(copEV))
	FE = AA(sort)(cop2lar(copFE))
	FV = AA(sort)([union(CAT([EV[e] for e in f])) for f in FE])

	ev_text = Dict{Vector{Int},Int}(collect(zip(EV, 1:length(EV))))
	fv_text = Dict{Vector{Int},Int}(collect(zip(FV, 1:length(FV))))

	outerspace, __others = separate_outer_atom(V, atoms)

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
export VIEWATOMS


# //////////////////////////////////////////////////////////////////////////////
function LAR_BATCHES(lar::Lar; properties::Properties=Properties())

	V  =lar.V
	EV=haskey(lar.C,:EV ) ? lar.C[:EV] : nothing
	FV=haskey(lar.C,:FV ) ? lar.C[:FV] : nothing

	properties["background_color"] = get(properties, "background_color", DEFAULT_LAR_BACKGROUND_COLOR)
	properties["line_color"] = get(properties, "line_color", DEFAULT_LINE_COLOR)
	properties["line_width"] = get(properties, "line_width", DEFAULT_LINE_WIDTH)

	properties["text_v_color"] = get(properties, "text_v_color", DEFAULT_TEXT_V_COLOR)
	properties["text_ev_color"] = get(properties, "text_ev_color", DEFAULT_TEXT_EV_COLOR)
	properties["text_fv_color"] = get(properties, "text_fv_color", DEFAULT_TEXT_FV_COLOR)

	# font sizes
	properties["v_fontsize" ] = get(properties, "v_fontsize", DEFAULT_V_FONTSIZE)
	properties["ev_fontsize"] = get(properties, "ev_fontsize", DEFAULT_EV_FONTSIZE)
	properties["fv_fontsize"] = get(properties, "fv_fontsize", DEFAULT_FV_FONTSIZE)

	Vtext  = haskey(lar.text,:V ) ? lar.text[:V ] : Dict( I => string(I) for I in eachindex(V ) )
	EVtext = haskey(lar.text,:EV) ? lar.text[:EV] : Dict( I => string(I) for I in eachindex(EV) )
	FVtext = haskey(lar.text,:FV) ? lar.text[:FV] : Dict( I => string(I) for I in eachindex(FV) )

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
				(I in used_vertices ? Vtext[I] : ""),
				center=ComputeCentroid([vertices[it] for it in [I]]),
				color=properties["text_v_color"],
				fontsize=properties["v_fontsize"]))
		end
	end

	# show edges
	if !isnothing(EV) && properties["text_ev_color"][4] > 0.0 && properties["ev_fontsize"] > 0
		for I in eachindex(EV)
			append!(batches, GLText(EVtext[I],
				center=ComputeCentroid([vertices[it] for it in EV[I]]),
				color=properties["text_ev_color"],
				fontsize=properties["ev_fontsize"]))
		end
	end

	# show faces
	if !isnothing(FV) && properties["text_fv_color"][4] > 0.0 && properties["fv_fontsize"] > 0
		for I in eachindex(FV)
			append!(batches,
				GLText(FVtext[I],
					center=ComputeCentroid([vertices[it] for it in FV[I]]),
					color=properties["text_fv_color"],
					fontsize=properties["fv_fontsize"]))
		end
	end

	return batches
end
export LAR_BATCHES

# //////////////////////////////////////////////////////////////////////////////
function VIEWEXPLODED(V, CVs::Vector{Cells}, FVs::Vector{Cells}, EVs::Vector{Cells}; scale_factor=1.2)

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
	# scrgiorgio: not sure this is right...
	if false
		begin
			v = []
			for (k,it) in enumerate(EXPLODECELLS(V, CVs, scale_factor=scale_factor))
				face_color = Point4d(Plasm.COLORS[(k-1)%12+1] - (rand(Float64, 4) * 0.1))
				face_color[4] = 1.0
				push!(v, PROPERTIES(it, Properties("line_width" => 3, "face_color" => face_color)))
			end
			VIEW(STRUCT(v))
		end
	end
end
export VIEWEXPLODED

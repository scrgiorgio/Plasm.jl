
# //////////////////////////////////////////////////////////////////////////////
function is_ray_intersecting_face(test_point::Vector{Float64}, face_points::Points)

	plane=create_plane(face_points)

	ray_origin = test_point
	ray_dir    = normalized([0, 0, 1.0]) # scrgiorgio: no need to normalize here... but just in case

	hit=plane_ray_intersection(ray_origin, ray_dir, plane)
	if isnothing(hit)
		return false
	end

	# i need to check if the hit is really inside the face
	# to do so I project all in 2d and use the 2d classify point
	projector   = project_points3d(face_points; double_check=true) # scrgiorgio: remove double check 
	hit         = projector(hit)
	face_points = projector(face_points)

	return classify_point(hit, BYROW(face_points), lar2cop(edges)) != "p_out"

end


# //////////////////////////////////////////////////////////////////////////////
function get_atom_internal_point(V::Points, FV::Cells, point::Vector{Float64}, normal::Vector{Float64})

	# note: I need to be sure I am around an internal/external position
	# I should move enough to avoid going in the error range (using 2-order of magnitute here)
	epsilon=LAR_DEFAULT_ERR*100

	p_test1 = point + epsilon * normal
	p_test2 = point - epsilon * normal

	is_internal1=(length([1 for face in FV if is_ray_intersecting_face(p_test1, V[:,face])]) % 2)==1
	is_internal2=(length([1 for face in FV if is_ray_intersecting_face(p_test2, V[:,face])]) % 2)==1

	if is_internal1 && !is_internal2
		return p_test1

	elseif is_internal2 && !is_internal1
		return p_test2

	end

	# ambiguous
	return nothing
end

# //////////////////////////////////////////////////////////////////////////////
function get_atom_internal_point(V::Points, FV::Cells) 

	# TODO: reintroduce space index to avoid O(N^2) complexity

	# here i am trying with all face centroids
	for face in FV
		face_points=V[:,face]
		ret=get_atom_internal_point(V,FV, get_centroid(face_points), plane_get_normal(create_plane(face_points)))
		if !isnothing(ret) return ret end
	end

	error("cannot find internal point")
end


# //////////////////////////////////////////////////////////////////////////////
function bool3d(assembly::Hpc, V::Points, EV::Cells, FE::Cells, CF::Cells)

	
	# separate outer atom (i.e. the one with the biggest diagonal)
	begin
		atoms,diags=[],[]
		for atom_faces in CF
			atom_ev=[[EV[E] for E in FE[F]] for F in atom_faces]
			atom_fv=[collect(Set(vcat([EV[E] for E in FE[F]]...))) for F in atom_faces]
			atom_fe=[collect(Set([E for E in FE[F]])) for F in atom_faces]
			push!(atoms, [atom_ev,atom_fv,atom_fe,atom_cf])

			# compute diagonal
			verts = sort(union(CAT(CAT(atom_ev)))) 
			bbox=collect([vec(it) for it in bbox_create(V[:, verts])])
			push!(diags, LinearAlgebra.norm(v2 - v1))

		end

		__value, outer_position = findmax(diags)
		outer      = filter(x -> x == atoms[outer_position], atoms)[1]
		atoms      = filter(x -> x != atoms[outer_position], atoms)
	end

	FV = [union(CAT([EV[e] for e in fe])) for fe in FE]

	# associate internal points to (original) faces of 3-cells
	begin
		lars = [LAR(it) for it in TOPOS(assembly)] 
		num_faces_per_atom = [LEN(it.C[:FV]) for it in lars] 
		cumulative = cumsum([0; num_faces_per_atom]) .+ 1
		fspans = collect(zip(cumulative[1:end-1], cumulative[2:end] .- 1))
		span(h) = [j for j = 1:length(fspans) if fspans[j][1] <= h <= fspans[j][2]]
	end

	# original assembly is needed for in/out ray test
	lar_assembly = LAR(assembly)
	assembly_V     = lar_assembly.V
	assembly_copEV = lar2cop(lar_assembly.C[:EV])
	assembly_copFV = lar2cop(lar_assembly.C[:FV])
	assembly_copFE = convert(ChainOp, assembly_copEV * assembly_copFV' .÷ 2)

	atoms_internal_points = []
	for (ev,fv,fe,cf) in atoms
		point, facenumber = get_atom_internal_point(V, fv, cf, copEV, copFE)
		push!(atoms_internal_points, point)
	end

	"""
	ret = {}
	ret[1, 1] = 1
	for (A, internal_point) in enumerate(atoms_internal_points) 
		faces=[]
		for assembly_face in assembly_FV
			if ray_intersect_face(internal_point, assembly_V, assembly_copEV, assembly_copFV, assembly_copFE, assembly_face)
				push!(faces, face)
			end
		end

		count = zeros(Int, length(fspans))
		for h in faces
			for s in span(h)
				count[s] += 1
			end
		end
		ret[A+1, 2:end] = count .% 2
	end
	
	return ret
	"""
end
export bool3d


# ///////////////////////////////////////////////////////////////
function subassembly3d(V_pass_through::Points, copEV::ChainOp, copFE::ChainOp, copCF::ChainOp, is_face_selected::Vector{Int})

	selected_faces = [face_index for face_index in 1:length(is_face_selected) if is_face_selected[face_index] ≠ 0]
	selected_edges = collect(Set([E for (F, E, Value) in zip(findnz(copFE)...) if is_face_selected[F] ≠ 0]))

	Fdict = Dict(zip(selected_faces, 1:length(selected_faces)))
	Edict = Dict(zip(selected_edges, 1:length(selected_edges)))
	
	triples_CF = hcat([[       C, Fdict[F], Value] for (C, F, Value) in zip(findnz(copCF)...) if is_face_selected[F] ≠ 0]...)
	triples_FE = hcat([[Fdict[F], Edict[E], Value] for (F, E, Value) in zip(findnz(copFE)...) if is_face_selected[F] ≠ 0]...)
	triples_EV = hcat([[Edict[E],        V, Value] for (E, V, Value) in zip(findnz(copEV)...) if     E in selected_edges]...)

	return (
		V_pass_through, 
		convert(SparseMatrixCSC{Int8,Int64}, sparse(triples_EV[1, :], triples_EV[2, :], triples_EV[3, :])),
		convert(SparseMatrixCSC{Int8,Int64}, sparse(triples_FE[1, :], triples_FE[2, :], triples_FE[3, :])), 
		convert(SparseMatrixCSC{Int8,Int64}, sparse(triples_CF[1, :], triples_CF[2, :], triples_CF[3, :]))
	)

end
export subassembly3d
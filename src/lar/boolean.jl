

# //////////////////////////////////////////////////////////////////////////////
"""
TODO: reintroduce space index to avoid O(N^2) complexity
"""
"""
function spaceindex_boolean(point3d::Array{Float64,1}, V::Points, CV::Cells)

	V=copy(V)
	CV=copy(CV)

	V = [V point3d]
	dim, idx = size(V)
	push!(CV, [idx, idx, idx])
	cellpoints = [V[:,cv]::Points for cv in CV]

	bboxes = [hcat(bbox_create(cell)...) for cell in cellpoints] 
	xboxdict = bbox_coord_intervals(1, bboxes)
	yboxdict = bbox_coord_intervals(2, bboxes)

	# xs,ys are IntervalTree type
	xs = IntervalTrees.IntervalMap{Float64,Array}()
	ys = IntervalTrees.IntervalMap{Float64,Array}()
	for (key, boxset) in xboxdict xs[tuple(key...)] = boxset end
	for (key, boxset) in yboxdict ys[tuple(key...)] = boxset end
	
	xcovers = bbox_covering(bboxes, 1, xs)
	ycovers = bbox_covering(bboxes, 2, ys)
	covers = [intersect(pair...) for pair in zip(xcovers, ycovers)]

	# remove each cell from its cover
	pointcover = setdiff(covers[end], [idx + 1])
	return pointcover[1:end-1]

end
"""

# //////////////////////////////////////////////////////////////////////////////
function get_ray_intersection_with_face(test_point::Vector{Float64}, face_points::Points)

	plane=create_plane(face_points)
	ray_origin = test_point
	ray_dir    = normalized([0, 0, 1.0]

	hit_3d=plane_ray_intersection(ray_origin, ray_dir), plane)

	if isnothing(hit_3d)
		return false
	end

	# i need to check if the hit is really inside the face
	# to do so I project all in 2d and use the 2d classify point

	center = compute_centroid(face_points)
  w=plane_get_normal(plane_create(face_points))
  v=normalized(cross(w,orthogonal_axis[argmin(w)]))
  u=normalized(cross(v,w))
	
	M=nothing# todo: apply  transation and rotation (!)

	hit_2d            = (M * hit_3d     )   [1:2,:]
	face_points_2d    = (M * face_points)   [1:2,:]
	
	classification=classify_point(hit_2d, BYROW(face_points_2d), lar2cop(edges))
	return classification!= "p_out"

end

# //////////////////////////////////////////////////////////////////////////////
function get_atom_internal_point(V::Points, FV::Cells, CF::Cells; epsilon=LAR_DEFAULT_ERR*100) 

	# for robustness, I am trying all atom faces. 
	# to be sure, I should get one internal and one external point
	for face in FV

		face_points=V[:,face]

		p0=get_centroid(face_points)

		# note the normal could be directed to the inside or the outside (==I don't want the faces to be coherently oriented)
		plane=create_plane(face_points)

		normal=plane_get_normal(plane)

		p_test1 = p0 + epsilon * normal
		p_test2 = p0 - epsilon * norma

		is_internal1=(length([for face in enumerate(CF) if get_ray_intersection_with_face(p_test1, face_points)]) % 2)==1
		is_internal2=(length([for face in enumerate(CF) if get_ray_intersection_with_face(p_test2, face_points)]) % 2)==1

		if is_internal1 && !is_internal2
			return p_test1

		elseif is_internal2 && !is_internal1
			return p_test2
		else
			# try with the next face
		end

	end

	error("cannot find internal point")
end



# //////////////////////////////////////////////////////////////////////////////
""" 
	V,copEV,copFE,copCF are coming from arrange3d
	input arguments are coming from arrage3d
"""

function bool3d(assembly::Hpc, V::Points, copEV::ChainOp, copFE::ChainOp, copCF::ChainOp)

	# sparse to dense 
	EV = cop2lar(copEV)
	FE = cop2lar(copFE) 
	CF = cop2lar(copCF)
	FV = [union(CAT([EV[e] for e in fe])) for fe in FE]

	# separate outer atom (i.e. the one with the biggest diagonal)
	begin
		atoms,diags=[],[]
		for k = 1:copCF.m
			cf=CF[k]
			ev=[[EV[e] for e in FE[f]] for f in cf]
			fv=[collect(Set(vcat([EV[e] for e in FE[f]]...))) for f in cf]
			fe=[collect(Set([e for e in FE[f]])) for f in cf]
			push!(atoms, [ev,fv,fe,cf])

			# compute diagonal
			verts = sort(union(CAT(CAT(ev)))) 
			bbox=collect([vec(it) for it in bbox_create(V[:, verts])])
			push!(diags, LinearAlgebra.norm(v2 - v1))

		end

		__value, outer_position = findmax(diags)
		outer      = filter(x -> x == atoms[outer_position], atoms)[1]
		atoms      = filter(x -> x != atoms[outer_position], atoms)
	end

	# VIEWATOMS(V,copEV,copFE,copCF, atoms; view_outer=true)

	# associate internal points to (original) faces of 3-cells
	lars = [LAR(it) for it in TOPOS(assembly)] 
	num_faces_per_atom = [LEN(it.C[:FV]) for it in lars] 
	cumulative = cumsum([0; num_faces_per_atom]) .+ 1
	fspans = collect(zip(cumulative[1:end-1], cumulative[2:end] .- 1))
	span(h) = [j for j = 1:length(fspans) if fspans[j][1] <= h <= fspans[j][2]]

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

	ret = {A:[]}
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
end
export bool3d



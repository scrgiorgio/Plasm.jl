
# //////////////////////////////////////////////////////////////////////////////
function get_atom_internal_point(Vcol::Points, FV::Cells, CF::Cells, copEV::ChainOp, copFE::ChainOp) 

	# first face
	f = CF[1]

	# first edge 
	e = findnz(copFE[f, :])[1][1] 

	# two faces incident on atom
	f1, f2 = findnz(copFE[f, :])[1][1:2] 

	# two verts incident
	v1, v2 = findnz(copEV[e, :])[1][1:2] 

	# enumerate atom faces
	fdict = Dict(zip(CF, 1:length(CF))) 
	vertices_f1 = FV[fdict[f1]]
	vertices_f2 = FV[fdict[f2]]

	# intersection should be an edge
	v1, v2 = intersect(vertices_f1, vertices_f2)

	# two triangles sharing a common edge 
	t1 = V_col[:, v1], V_col[:, v2], V_col[:, [v for v in vertices_f1 if v ≠ v1 && v ≠ v2][1]]
	t2 = V_col[:, v2], V_col[:, v1], V_col[:, [v for v in vertices_f2 if v ≠ v1 && v ≠ v2][1]]
	n1 = LinearAlgebra.normalize(cross(t1[2] - t1[1], t1[3] - t1[1]))
	n2 = LinearAlgebra.normalize(cross(t2[2] - t2[1], t2[3] - t2[1]))

	# mean point
	p0 = (V_col[:, v1] + V_col[:, v2]) ./ 2 
	
	# mean normal, I am assuming that moving a little from the border following the normal I will get an internal point
	#   which could be wrong for very thin solids
	n = (n1 + n2)/2.0 
	above = p0 + ϵ * n
	below = p0 - ϵ * n

	k1, k2 = 0, 0
	for (f, face) in enumerate(CF)
		hit1 = get_ray_intesection_with_plane(above, V_col[:, FV[f]])
		hit2 = get_ray_intesection_with_plane(below, V_col[:, FV[f]])
		if point_in_face_3d(hit1, V_col, copEV, copFE, face) k1 += 1 end 
		if point_in_face_3d(hit2, V_col, copEV, copFE, face) k2 += 1 end
	end

	if (k1 % 2) == 1
		return above

	elseif (k2 % 2) == 1
		return below

	else
		error("tertium non datur")
	end
end


# //////////////////////////////////////////////////////////////////
function point_in_face_3d(hit, V_col::Points, copEV::ChainOp, copFE::ChainOp, face::Int)

	if isnothing(hit)
		return false
	end

	fv, edges = find_vcycle(copEV, copFE, face)

	V_col = V_col[:, fv]
	
	first=V_col[:, 1]

	hit   = hit .- first
	V_col     = V_col .- first
	u, v = edges[1]
	__z, w = [[z, w] for (z, w) in edges if z == v][1]

	v1 = V_col[:, u] - V_col[:, v]
	v2 = V_col[:, w] - V_col[:, v]
	v3 = cross(v2, v1)
	
	M = [v1 v2 v3]
	projected = inv(M) * [hit point]

	return classify_point(projected[1:2, end], BYROW(projected[1:2, 1:end-1]), lar2cop(edges))()!= "p_out"

end



# //////////////////////////////////////////////////////////////////////////////
function spaceindex_boolean(point3d::Array{Float64,1}, V_col::Points, CV::Cells)

	V_col=copy(V_col)
	CV=copy(CV)

	V_col = [V_col point3d]
	dim, idx = size(V_col)
	push!(CV, [idx, idx, idx])
	cellpoints = [V_col[:,cv]::Points for cv in CV]

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

# //////////////////////////////////////////////////////////////////////////////
function get_ray_intesection_with_plane(point3d, face_points::Points)

	p0=face_points[:, 1]
	p1=face_points[:, 2]
	p2=face_points[:, 3]

	v1=p1-p0
	v2=p2-p0
	n = LinearAlgebra.normalize(cross(v1, v2))

	l= [0, 0, 1.0]
	denom = dot(n, l)
	if (abs(denom) > LAR_DEFAULT_ERR)
		t = dot(p0 - point3d, n) / denom
		if t > 0
			return point3d + t * l   
		end
	end

	# no intersection or behind
	return nothing
end

# //////////////////////////////////////////////////////////////////////////////
function get_faces_intersecting_ray(point3d, V_col::Points, copEV::ChainOp, coFV::ChainOp, copFE::ChainOp)
	ret = Int64[]
	candidates=spaceindex_boolean(point3d, V_col, FV)
	for face in candidates
		hit = get_ray_intesection_with_plane(point3d, V_col[:, FV[face]])
		if point_in_face_3d(hit, V_col, copEV, copFE, face)
			push!(ret, face)
		end
	end
	return ret 
end

# //////////////////////////////////////////////////////////////////////////////
function bool3d(V_col::Points, copEV::ChainOp, copFE::ChainOp, copCF::ChainOp, assembly::Hpc)

	lar = LAR(assembly)
	assembly_V     = lar.V
	assembly_copEV = lar2cop(lar.C[:EV])
	assembly_copFV = lar2cop(lar.C[:FV])
	assembly_copFE = convert(ChainOp, assembly_copEV * assembly_copFV' .÷ 2)

	# atoms
	EV = cop2lar(copEV)
	FE = cop2lar(copFE) 
	CF = cop2lar(copCF)
	FV = [union(CAT([EV[e] for e in fe])) for fe in FE]

	atoms,diags=[],[]
	for k = 1:copCF.m
		cf=CF[k]
		ev=[[EV[e] for e in FE[f]] for f in cf]
		fv=[collect(Set(vcat([EV[e] for e in FE[f]]...))) for f in cf]
		fe=[collect(Set([e for e in FE[f]])) for f in cf]
		push!(atoms, [ev,fv,fe,cf])

		# compute diagonal
		verts = sort(union(CAT(CAT(ev)))) 
		bbox=collect([vec(it) for it in bbox_create(V_col[:, verts])])
		push!(diags, LinearAlgebra.norm(v2 - v1))

	end

	__value, outer_position = findmax(diags)
	outer      = filter(x -> x == atoms[outer_position], atoms)[1]
	others     = filter(x -> x != atoms[outer_position], atoms)

	# VIEWATOMS(V,copEV,copFE,copCF, atoms; view_outer=true)

	atoms_inner_points = []
	for (ev,fv,fe,cf) in others
		point, facenumber = get_atom_internal_point(Vcol, fv, cf, copEV, copFE)
		push!(atoms_inner_points, point)
	end

	# associate internal points to (original) faces of 3-cells
	lars = [LAR(it) for it in TOPOS(assembly)] 

	num_faces_per_atom = [LEN(it.C[:FV]) for it in lars] 
	cumulative = cumsum([0; num_faces_per_atom]) .+ 1
	fspans = collect(zip(cumulative[1:end-1], cumulative[2:end] .- 1))
	span(h) = [j for j = 1:length(fspans) if fspans[j][1] <= h <= fspans[j][2]]

	ret = BitArray(undef, length(atoms_inner_points) + 1, length(fspans) + 1)
	ret[1, 1] = 1
	for (A, atom_inner_point) in enumerate(atoms_inner_points) 
		faces = get_faces_intersecting_ray(atom_inner_point, assembly_V, assembly_copEV, assembly_copFV, assembly_copFE)
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



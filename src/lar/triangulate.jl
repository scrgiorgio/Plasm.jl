

# //////////////////////////////////////////////////////////////////////
function constrained_triangulation2D(V::Points, EV::Cells)
	triin = Triangulate.TriangulateIO()
	triin.pointlist = V # scrgiorgio: this is by-col representation as LAR
	triin.segmentlist = hcat(EV...)
	(triout, __vorout) = Triangulate.triangulate("pQ", triin)  # exec triangulation
	return Array{Int64,1}[c[:] for c in eachcol(triout.trianglelist)]
end

# //////////////////////////////////////////////////////////////////////
function TRIANGULATE2D(V::Points, EV::Cells)
	num_vertices=size(V,2)
	copEV = lar2cop(EV)
	V_row = BYROW(V)
	trias = constrained_triangulation2D(V, EV)
	ret = Array{Int,1}[]
	for (u, v, w) in trias
		centroid = (V_row[u, :] + V_row[v, :] + V_row[w, :]) ./ 3
		if point_in_face(centroid, V_row, copEV)
			push!(ret, [u, v, w])
		end
	end
	return ret
end
export TRIANGULATE2D

# //////////////////////////////////////////////////////////////////////
""" input old LAR consistent data; output triangulated_faces """
function LAR2TRIANGLES(V::Points, EV::Cells, FV::Cells, FE::Cells;err = LAR_DEFAULT_ERR)

	triangles_per_face = Vector{Any}(undef, length(FE))

	for edges_idxs in FE
		edge_num = length(edges_idxs)
		fv, edges = find_cycle2(EV, FE, f)
		# look for independent vector triple
		points = V[:, fv]
		vmap = Dict(zip(fv, 1:length(fv))) # vertex map
		mapv = Dict(zip(1:length(fv), fv)) # inverse vertex map
		edges = [[vmap[A], vmap[B]] for (A, B) in edges]
		v1 = LinearAlgebra.normalize(points[2, :] - points[1, :])
		v2 = [0, 0, 0]
		v3 = [0, 0, 0]
		i = 3
		while -err < LinearAlgebra.norm(v3) < err
			v2 = LinearAlgebra.normalize(points[i, :] - points[1, :])
			v3 = LinearAlgebra.cross(v1, v2)
			i = i % size(points, 1) + 1
		end
		
		# independent vector triple in face f 
		M = [v1 v2 v3]
		projected = BYCOL((points*M)[:, 1:2])
		triangles_per_face[f] = [[mapv[v] for v in t] for t in constrained_triangulation2D(projected, edges)]
	end

	return triangles_per_face
end
export LAR2TRIANGLES

# //////////////////////////////////////////////////////////////////////
""" From  topology to cells (1D chains, 2D chains, breps of 3D chains) """
function pols2tria(V::Points, copEV::ChainOp, copFE::ChainOp, copCF) 

	# historically I need to access vertices by row
	V_row = BYROW(V)

	triangles_per_face = Vector{Any}(undef, copFE.m)

	for f in 1:copFE.m
		edges_idxs = copFE[f, :].nzind
		edge_num = length(edges_idxs)
		edges = zeros(Int, edge_num, 2)
		fv, edges = find_vcycle_v1(copEV, copFE, f) # scrgiorgio: find_vcycle_v2 DOES NOT WORK! (why???)
		if fv ≠ []
			vs = V_row[fv, :]
			v1 = LinearAlgebra.normalize(vs[2, :] - vs[1, :])
			v2 = [0, 0, 0]
			v3 = [0, 0, 0]
			i = 3
			while LinearAlgebra.norm(v3) < abs(LAR_DEFAULT_ERR)
				v2 = LinearAlgebra.normalize(vs[i, :] - vs[1, :])
				v3 = LinearAlgebra.cross(v1, v2)
				i = i % size(vs, 1) + 1
			end
			M = reshape([v1; v2; v3], 3, 3)
			vs = (vs*M)[:, 1:2]
			v = convert(Points, vs'[1:2, :])
			vmap = Dict(zip(fv, 1:length(fv))) # vertex map
			mapv = Dict(zip(1:length(fv), fv)) # inverse vertex map
			triangles_per_face[f] = [[mapv[v] for v in tria] for tria in TRIANGULATE2D(v, edges)]
		end
	end
	triangles_per_face = convert(Vector{Cells}, triangles_per_face)

	# polygonal face fragments
	EVs = FV2EVs(copEV, copFE) 
	FVs = convert(Array{Cells}, triangles_per_face)
	CVs = []
	for cell in 1:copCF.m
		obj = []
		for f in copCF[cell, :].nzind
			append!(obj, triangles_per_face[f])
		end
		push!(CVs, obj)
	end
	return V, convert(Vector{Cells}, CVs), FVs, EVs
end
export pols2tria

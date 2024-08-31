# //////////////////////////////////////////////////////////////////////
function constrained_triangulation2D(V::Points, EV::Cells)
	triin = Triangulate.TriangulateIO()
	triin.pointlist = V # scrgiorgio: by-col representation as LAR
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
function LAR2TRIANGLES(V::Points, EV::Cells, FV::Cells, FE::Cells;err = 1e-8)

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

	triangles_per_face = Vector{Any}(undef, length(FE))

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
function pols2tria(V, copEV, copFE, copCF) # W by columns
	V_row = BYROW(V)

	triangles_per_face = Vector{Any}(undef, copFE.m)

	for f in 1:copFE.m
		if f % 10 == 0
			print(".")
		end
		edges_idxs = copFE[f, :].nzind
		edge_num = length(edges_idxs)
		edges = zeros(Int, edge_num, 2)
		fv, edges = vcycle(copEV, copFE, f)
		if fv ≠ []
			vs = V_row[fv, :]
			v1 = LinearAlgebra.normalize(vs[2, :] - vs[1, :])
			v2 = [0, 0, 0]
			v3 = [0, 0, 0]
			err = 1e-8
			i = 3
			while -err < LinearAlgebra.norm(v3) < err
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
	return V, CVs, FVs, EVs
end
export pols2tria

# ///////////////////////////////////////////////////////////////
## Fs is the signed coord vector of a subassembly
## the logic is to compute the corresponding reduced coboundary matrices
## and finally call the standard method of the function.
function pols2tria(V, copEV, copFE, copCF, Fs)
	# make copies of coboundary operators

	# compute the reduced copCF
	CFtriples = findnz(copCF)
	triples = [triple for triple in zip(CFtriples...)]
	newtriples = [(row, col, val) for (row, col, val) in triples if Fs[col] ≠ 0]
	newF = [k for (k, f) in enumerate(Fs) if Fs[k] ≠ 0]
	fdict = Dict(zip(newF, 1:length(newF)))
	triples = hcat([[row, fdict[col], val] for (row, col, val) in newtriples]...)
	newCF = sparse(triples[1, :], triples[2, :], triples[3, :])
	copCF = convert(SparseMatrixCSC{Int8,Int64}, newCF)

	# compute the reduced copFE
	FEtriples = findnz(copFE)
	triples = [triple for triple in zip(FEtriples...)]
	newtriples = [(row, col, val) for (row, col, val) in triples if Fs[row] ≠ 0]
	newF = [k for (k, f) in enumerate(Fs) if Fs[k] ≠ 0]
	newcol = collect(Set([col for (row, col, val) in newtriples]))
	facedict = Dict(zip(newF, 1:length(newF)))
	edgedict = Dict(zip(newcol, 1:length(newcol)))
	triples = hcat([[facedict[row], edgedict[col], val] for (row, col, val) in newtriples]...)
	newFE = sparse(triples[1, :], triples[2, :], triples[3, :])
	copFE = convert(SparseMatrixCSC{Int8,Int64}, newFE)

	# compute the reduced copEV
	EVtriples = findnz(copEV)
	triples = [triple for triple in zip(EVtriples...)]
	newtriples = [(row, col, val) for (row, col, val) in triples if row in keys(edgedict)]
	triples = hcat([[edgedict[row], col, val] for (row, col, val) in newtriples]...)
	newEV = sparse(triples[1, :], triples[2, :], triples[3, :])
	copEV = convert(SparseMatrixCSC{Int8,Int64}, newEV)

	# finally compute the cells, faces, and edges of subassembly
	return pols2tria(V, copEV, copFE, copCF)
end
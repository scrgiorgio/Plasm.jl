



# //////////////////////////////////////////////////////////////////////////////
function tgw_containment_test(triangle, f, V, FV, EV, FE)::Bool
	""" Triangle containment test of checkpoint into `f`; used in TGW algorithm """
	# affine mapping computation to z=0 plane
	v1, v2, v3 = triangle
	u = V[:, v2] - V[:, v1]
	v = V[:, v3] - V[:, v1]
	w = cross(u, v)
	T = [1.0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
	T[1:3, 4] = -V[:, v1]
	R = [1.0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
	R[1:3, 1:3] = [u v w]
	R = R'
	mapping = R * T

	trianglepoints = [[V[:, v1] V[:, v2] V[:, v3]]; ones(3)']
	pts = mapping * trianglepoints
	points2D = [pts[r, :] for r = 1:size(pts, 1)
							if !(all(pts[r, :] .== 0) || all(pts[r, :] .== 1))]
	p2D = hcat(points2D...)'
	checkpoint = 0.4995 .* p2D[:, 1] + 0.4995 .* p2D[:, 2] + 0.001 .* p2D[:, 3]

	cellpoints = [V[:, FV[f]]; ones(length(FV[f]))']
	points = mapping * cellpoints
	verts2D = [points[r, :] for r = 1:size(points, 1)
						 if !(all(points[r, :] .== 0) || all(points[r, :] .== 1))]
	P2D = hcat(verts2D...)

	vdict = Dict(collect(zip(FV[f], 1:length(FV[f]))))
	celledges = [[vdict[v] for v in EV[e]] for e in FE[f]]
	inner = point_in_face(checkpoint, P2D::Points, lar2cop(celledges)::ChainOp)
end


# //////////////////////////////////////////////////////////////////////////////
function tgw_triangle_ordering(triangles, V)
	""" Correct ordering of triangles (hence faces) about each boundary edge """
	normals = []
	v1, v2, v3 = triangles[1]
	if v1 > v2
		v1, v2 = v2, v1
	end
	e3 = LinearAlgebra.normalize(V[:, v2] - V[:, v1])
	e1 = LinearAlgebra.normalize(V[:, v3] - V[:, v1])
	e2 = LinearAlgebra.normalize(cross(e1, e3))
	basis = [e1 e2 e3]
	transform = inv(basis)

	angles = []
	for (v1, v2, v3) in triangles
		w1 = LinearAlgebra.normalize(V[:, v3] - V[:, v1])
		w2 = transform * w1
		w3 = cross([0, 0, 1], w2)
		push!(normals, w3)
	end
	for k = 1:length(normals)
		angle = atan(normals[k][2], normals[k][1])
		push!(angles, angle)
	end
	pairs = sort(collect(zip(angles, 1:length(triangles))))
	order = [k for (angle, k) in pairs]
	return order
end


# //////////////////////////////////////////////////////////////////////////////
function tgw_compute_ordered_fan(hinge::Int, bd1, V::Points, FV::Cells, EV::Cells, FE::Cells)

	""" Component of TGW in 3D (Pao);  return an ordered `fan` of 2-cells """
	cells = SparseArrays.findnz(bd1)[1]
	triangles = []

	function area(v1, v2, v3)
		u = V[:, v2] - V[:, v1]
		v = V[:, v3] - V[:, v1]
		out = LinearAlgebra.norm(cross(u, v)) # actually, to be divided by two
		return out
	end

	for f in cells
		v1, v2 = EV[hinge]
		index = findfirst(v -> (area(v1, v2, v) ≠ 0), FV[f])
		v3 = FV[f][index]

		# test if [v1,v2,v3] interior to f
		while true
			if tgw_containment_test([v1, v2, v3], f, V, FV, EV, FE)
				push!(triangles, [v1, v2, v3])
				break
			else
				index = findnext(v -> (area(v1, v2, v) ≠ 0), FV[f], index + 1)
				v3 = FV[f][index]
			end
		end
	end
	order = tgw_triangle_ordering(triangles, V)
	return [cells[index] for index in order]
end


# //////////////////////////////////////////////////////////////////////////////
function tgw_next(cycle, pivot, marks)
	
	""" Utility function for TGW in 3D """
	len = length(cycle)
	ind = findfirst(x -> x == pivot, cycle)[1]
	nextIndex = ind == len ? 1 : ind + 1
	#if marks[nextIndex]==2
	#  nextIndex = ind==1 ? len : ind-1
	#end
	return cycle[nextIndex][1]
end

# //////////////////////////////////////////////////////////////////////////////
function tgw_prev(cycle, pivot, marks)
	""" Utility function for TGW in 3D """
	len = length(cycle)
	ind = findfirst(x -> x == pivot, cycle)[1] #  ind is address of pivot in cycle
	nextIndex = ind == 1 ? len : ind - 1
	return cycle[nextIndex][1]
end

# //////////////////////////////////////////////////////////////////////////////
function tgw_build_copCF(V_row, rcopEV, rcopFE)
	
	""" TGW algorithm implementation (pao) """

	# G&F -> Pao data structures
	V = BYCOL(V_row)
	num_vertices=size(V, 2)
	EV = cop2lar(rcopEV)
	fe = cop2lar(rcopFE)
	fv = [union([EV[e] for e in fe[f]]...) for f = 1:length(fe)]
	FV = convert(Cells, fv)
	copFE = rcopFE    # alias useful for algorithm testing 
	VV = [[index] for index = 1:num_vertices]

	copEF = transpose(copFE)
	FE = [SparseArrays.findnz(copFE[k, :])[1] for k = 1:size(copFE, 1)]
	# Initializations
	m, n = size(copEF)
	marks = zeros(Int, n)
	I = Int[]
	J = Int[]
	W = Int[]
	jcol = 0
	choose(marks) = findall(x -> x < 2, marks)[end]

	# Main loop (adding one copFC's column stepwise)
	# while sum(marks) < 2n # no robust condition ... make better
	while (!all(value == 2 for value in marks) && (sum(marks) <= 2n))
		# to don't loop
		# select a (d−1)-cell, "seed" of the column extraction
		σ = choose(marks)
		if marks[σ] == 0
			cd1 = sparsevec([σ], Int[1], n)
		elseif marks[σ] == 1
			cd1 = sparsevec([σ], Int[-1], n)
		end
		# compute boundary cd2 of seed cell
		cd2 = copEF * cd1
		# loop until (boundary) cd2 becomes empty
		while nnz(cd2) ≠ 0
			corolla = sparsevec([], Int[], m)
			# for each “hinge” τ cell
			for τ ∈ (.*)(SparseArrays.findnz(cd2)...)
				#compute the  coboundary
				tau = sparsevec([abs(τ)], Int[sign(τ)], m)  # ERROR: index out of bound here!
				bd1 = transpose(transpose(tau) * copEF)
				cells2D = SparseArrays.findnz(bd1)[1]
				# compute the  support
				inters = intersect(cells2D, SparseArrays.findnz(cd1)[1])
				if inters ≠ []
					pivot = inters[1]
				else
					error("no pivot")
				end
				# compute the new adj cell
				fan = tgw_compute_ordered_fan(abs(τ), bd1, V, FV, EV, FE) 
				if τ > 0
					adj = tgw_next(fan, pivot, marks)
				elseif τ < 0
					adj = tgw_prev(fan, pivot, marks)
				end
				# orient adj
				if copEF[abs(τ), adj] ≠ copEF[abs(τ), pivot]
					corolla[adj] = cd1[pivot]
				else
					corolla[adj] = -(cd1[pivot])
				end
			end
			# insert corolla cells in current cd1
			for (k, val) in zip(SparseArrays.findnz(corolla)...)
				cd1[k] = val
			end
			# compute again the boundary of cd1
			cd2 = copEF * cd1
		end
		for σ ∈ SparseArrays.findnz(cd1)[1]
			# update the counters of used cells
			marks[σ] += 1
		end
		# append a new column to [∂d+]
		# copFC += cd1
		rows, vals = SparseArrays.findnz(cd1)
		jcol += 1
		append!(I, rows)
		append!(J, [jcol for k = 1:nnz(cd1)])
		append!(W, vals)
	end
	copCF = sparse(J, I, W)
	# @show copCF;
	return copCF
end


#//////////////////////////////////////////////////////////////////////////////
# REMARK: will not work when vs[:,1:3] are aligned !!!!  TODO: fix 
function fragment_submanifold_mapping(vs)
	""" Compute the map from vs (at least three) to z=0 """
	u1 = vs[2, :] - vs[1, :]
	u2 = vs[3, :] - vs[1, :]
	u3 = LinearAlgebra.cross(u1, u2)
	T = Matrix{Float64}(LinearAlgebra.I, 4, 4)
	T[4, 1:3] = -vs[1, :]
	M = Matrix{Float64}(LinearAlgebra.I, 4, 4)
	M[1:3, 1:3] = [u1 u2 u3]
	return T * M
end


#//////////////////////////////////////////////////////////////////////////////
function fragment_face_intersection(V::Points, EV::ChainOp, face::Chain; err=LAR_DEFAULT_ERR)
	vs = buildFV(EV, face)     # EV::ChainOp, face::Chain
	retV = Points(undef, 0, 3)
	visited_verts = []
	for i in 1:length(vs)
		o = V[vs[i], :]
		j = i < length(vs) ? i + 1 : 1
		d = V[vs[j], :] - o    # vertex in local coordinates
		if !(-err < d[3] < err)
			alpha = -o[3] / d[3]
			if -err <= alpha <= 1 + err
				p = o + alpha * d
				if -err < alpha < err || 1 - err < alpha < 1 + err
					if !(is_visited_vertex(p, visited_verts))
						push!(visited_verts, p)
						retV = [retV; reshape(p, 1, 3)]
					end
				else
					retV = [retV; reshape(p, 1, 3)]
				end
			end
		end
	end
	vnum = size(retV, 1)
	if vnum == 1
		vnum = 0
		retV = Points(undef, 0, 3)
	end
	enum = (÷)(vnum, 2)
	retEV = spzeros(Int8, enum, vnum)
	for i in 1:enum
		retEV[i, 2*i-1:2*i] = [-1, 1]
	end
	retV, retEV
end



# //////////////////////////////////////////////////////////////////////////////
function fragment_merge_vertices(V::Points, EV::ChainOp, FE::ChainOp; err=LAR_DEFAULT_ERR)
	""" Task to iteratively add new local components to the global 2-skeleton """
	vertsnum = size(V, 1)
	edgenum = size(EV, 1)
	facenum = size(FE, 1)
	newverts = zeros(Int, vertsnum)
	# KDTree constructor needs an explicit array of Float64
	V = Matrix(V)
	kdtree = KDTree(BYROW(V))
	# remove vertices congruent to a single representative
	todelete = []
	i = 1
	for vi in 1:vertsnum
		if !(vi in todelete)
			nearvs = inrange(kdtree, V[vi, :], err)
			newverts[nearvs] .= i
			nearvs = setdiff(nearvs, vi)
			todelete = union(todelete, nearvs)
			i = i + 1
		end
	end
	nV = V[setdiff(collect(1:vertsnum), todelete), :]
	# translate edges to take congruence into account
	edges = Array{Tuple{Int,Int},1}(undef, edgenum)
	oedges = Array{Tuple{Int,Int},1}(undef, edgenum)
	for ei in 1:edgenum
		v1, v2 = EV[ei, :].nzind
		edges[ei] = Tuple{Int,Int}(sort([newverts[v1], newverts[v2]]))
		oedges[ei] = Tuple{Int,Int}(sort([v1, v2]))
	end
	nedges = union(edges)
	# remove edges of zero length
	nedges = filter(t -> t[1] != t[2], nedges)
	nedgenum = length(nedges)
	nEV = spzeros(Int8, nedgenum, size(nV, 1))
	etuple2idx = Dict{Tuple{Int,Int},Int}()
	for ei in 1:nedgenum
		begin
			nEV[ei, collect(nedges[ei])] .= 1
			nEV
		end
		etuple2idx[nedges[ei]] = ei
	end
	for e in 1:nedgenum
		v1, v2 = findnz(nEV[e, :])[1]
		nEV[e, v1] = -1
		nEV[e, v2] = 1
	end
	# compute new faces to take congruence into account
	faces = [[
		map(x -> newverts[x], FE[fi, ei] > 0 ? oedges[ei] : reverse(oedges[ei]))
		for ei in FE[fi, :].nzind
	] for fi in 1:facenum]

	visited = []
	function filter_fn(face)
		verts = []
		map(e -> verts = union(verts, collect(e)), face)
		verts = Set(verts)
		if !(verts in visited)
			push!(visited, verts)
			return true
		end
		return false
	end
	nfaces = filter(filter_fn, faces)
	nfacenum = length(nfaces)
	nFE = spzeros(Int8, nfacenum, size(nEV, 1))
	for fi in 1:nfacenum
		for edge in nfaces[fi]
			ei = etuple2idx[Tuple{Int,Int}(sort(collect(edge)))]
			nFE[fi, ei] = sign(edge[2] - edge[1])
		end
	end
	return Points(nV), nEV, nFE
end

# //////////////////////////////////////////////////////////////////////////////
function fragment_single_face(V, EV::ChainOp, FE::ChainOp, sp_idx, sigma)
	vs_num = size(V, 1)
	# 2D transformation of `sigma` face
	sigmavs = (abs.(FE[sigma:sigma, :])*abs.(EV))[1, :].nzind # sigma vertex indices
	sV = V[sigmavs, :]
	sEV = EV[FE[sigma, :].nzind, sigmavs]
	M = fragment_submanifold_mapping(sV)
	tV = ([V ones(vs_num)]*M)[:, 1:3]  # folle convertire *tutti* i vertici
	sV = tV[sigmavs, :]
	# `sigma` face intersection with faces in `sp_idx[sigma]`, i.e., in `bigpi`
	for i in sp_idx[sigma]
		tmpV, tmpEV = fragment_face_intersection(tV, EV, FE[i, :]) # va in loop qui dentro
		sV, sEV = skel_merge(sV, sEV, tmpV, tmpEV)
	end
	# computation of 2D arrangement of sigma face
	sV = sV[:, 1:2]
	nV, nEV, nFE = planar_arrangement(BYCOL(sV), sEV, sparsevec(ones(Int8, length(sigmavs))))
	nV = BYROW(nV)
	nvsize = size(nV, 1)
	# return each 2D complex in 3D
	nV = [nV zeros(nvsize) ones(nvsize)] * inv(M)[:, 1:3]
	return nV, nEV, nFE
end


# //////////////////////////////////////////////////////////////////////////////
function fragment_compute_copFE(V::Points, copFV::ChainOp, copEV::ChainOp; convex=true::Bool, exterior=false::Bool)::ChainOp

	temp = copFV * copEV'
	I, J, Val = Int64[], Int64[], Int8[]
	for j = 1:size(temp, 2)
		for i = 1:size(temp, 1)
			if temp[i, j] == 2
				push!(I, i)
				push!(J, j)
				push!(Val, 1)
			end
		end
	end
	copFE = SparseArrays.sparse(I, J, Val)
	if !convex
		copFE = fix_redundancy(copFE, copFV, copEV)
	end

	EV = [findnz(copEV[k, :])[1] for k = 1:size(copEV, 1)]
	copEV = sparse(cop_coboundary_0(EV))
	for f = 1:size(copFE, 1)
		chain = findnz(copFE[f, :])[1]#	dense
		cycle = spzeros(Int8, copFE.n)#	sparse

		edge = findnz(copFE[f, :])[1][1]
		sign = 1
		cycle[edge] = sign
		chain = setdiff(chain, edge)
		while chain != []
			boundary = sparse(cycle') * copEV
			_, vs, vals = findnz(dropzeros(boundary))

			rindex = vals[1] == 1 ? vf = vs[1] : vf = vs[2]
			r_boundary = spzeros(Int8, copEV.n)#	sparse
			r_boundary[rindex] = 1
			r_coboundary = copEV * r_boundary
			r_edge = intersect(findnz(r_coboundary)[1], chain)[1]
			r_coboundary = spzeros(Int8, copEV.m)#	sparse
			r_coboundary[r_edge] = EV[r_edge][1] < EV[r_edge][2] ? 1 : -1

			lindex = vals[1] == -1 ? vi = vs[1] : vi = vs[2]
			l_boundary = spzeros(Int8, copEV.n)#	sparse
			l_boundary[lindex] = -1
			l_coboundary = copEV * l_boundary
			l_edge = intersect(findnz(l_coboundary)[1], chain)[1]
			l_coboundary = spzeros(Int8, copEV.m)#	sparse
			l_coboundary[l_edge] = EV[l_edge][1] < EV[l_edge][2] ? -1 : 1

			if r_coboundary != -l_coboundary  # false iff last edge
				# add edge to cycle from both sides
				rsign = rindex == EV[r_edge][1] ? 1 : -1
				lsign = lindex == EV[l_edge][2] ? -1 : 1
				cycle = cycle + rsign * r_coboundary + lsign * l_coboundary
			else
				# add last (odd) edge to cycle
				rsign = rindex == EV[r_edge][1] ? 1 : -1
				cycle = cycle + rsign * r_coboundary
			end
			chain = setdiff(chain, findnz(cycle)[1])
		end
		for e in findnz(copFE[f, :])[1]
			copFE[f, e] = cycle[e]
		end
	end

	pdim=size(V, 1)

	if exterior &&  pdim == 2
		# put matrix in form: first row outer cell; with opposite sign )
		V_row  = BYROW(V)

		outer = get_external_cycle(V_row, copEV, copFE)
		copFE = [-copFE[outer:outer, :]; copFE[1:outer-1, :]; copFE[outer+1:end, :]]
		# induce coherent orientation of matrix rows (see examples/orient2d.jl)
		for k = 1:size(copFE, 2)
			spcolumn = findnz(copFE[:, k])
			if sum(spcolumn[2]) != 0
				row = spcolumn[1][2]
				sign = spcolumn[2][2]
				copFE[row, :] = -sign * copFE[row, :]
			end
		end
		return copFE
	else
		return copFE
	end
end


# ///////////////////////////////////////////////////////////
function fragment_all_faces(lar::Lar)

	V=lar.V
	EV=lar.C[:EV]
	FV=lar.C[:FV]

	# scrgiorgio: if I use here simply `lar2cop` it does not work
	# the cop_XXX function do some magic with orientation
	if true
		copEV = cop_coboundary_0(EV)
		copFE = fragment_compute_copFE(V, lar2cop(FV), lar2cop(EV), convex=true, exterior=false)
	else
		copEV = lar2cop(EV)
		copFV = lar2cop(FV)
		copFE = (copFV * copEV') .÷ Int8(2)
	end
	
	# historically arrangement works internally by using by-row vertices
	V_row = BYROW(V)
	fs_num = size(copFE, 1)
	# strange but necessary cycle of computations to get FV::Cells algebraically
	FV = (abs.(copFE) * abs.(copEV)) .÷ 2
	FV = convert(ChainOp, FV)
	sp_idx = spaceindex(V_row, cop2lar(FV))
	rV = Points(undef, 0, 3)
	rEV = SparseArrays.spzeros(Int8, 0, 0)
	rFE = SparseArrays.spzeros(Int8, 0, 0)
	depot_V = Array{Array{Float64,2},1}(undef, fs_num)
	depot_EV = Array{ChainOp,1}(undef, fs_num)
	depot_FE = Array{ChainOp,1}(undef, fs_num)
	for sigma in 1:fs_num
		print(sigma, "/", fs_num, "\r")
		nV, nEV, nFE = fragment_single_face(V_row, copEV, copFE, sp_idx, sigma)
		depot_V[sigma] = nV
		depot_EV[sigma] = nEV
		depot_FE[sigma] = nFE
	end
	rV = vcat(depot_V...)
	rEV = SparseArrays.blockdiag(depot_EV...)
	rFE = SparseArrays.blockdiag(depot_FE...)
	rV, rcopEV, rcopFE = fragment_merge_vertices(rV, rEV, rFE)
	return rV, rcopEV, rcopFE
end

# //////////////////////////////////////////////////////////////////////////////

function ARRANGE3D(lar::Lar; debug_mode=false)::Lar
	""" Main function of arrangement pipeline """
	rV, rcopEV, rcopFE=fragment_all_faces(lar)

  EV = cop2lar(rcopEV)
  FE = cop2lar(rcopFE) 
  FV = [union(CAT([EV[E] for E in fe])) for fe in FE]

  ret = Lar(BYCOL(rV),Dict(
		:EV => EV, 
		:FE => FE, 
		:FV => FV
	))

	if debug_mode
		VIEWCOMPLEX(ret, explode=[1.0,1.0,1.0], show=["V","EV","Vtext"], title="arrange3d / 3d ALL faces")
	end  

	if LAR_ARRANGE3D_USE_EXPERIMENTAL
		ret=SIMPLIFY(ret) # not sure this is needed
		ret, cycles=explode_cycles(ret)
		ret.C[:FV]=compute_FV(ret)
		ret.C[:CF]=lar_find_atoms(ret.V, cycles, debug_mode=debug_mode)

		# ret.C[:CV]=compute_CV(ret,is_convex=true) dont think this is needed
		CHECK(ret)
	else
		# broken, fails to find atoms in case of disconnected components
		rcopCF = tgw_build_copCF(rV, rcopEV, rcopFE)
		ret.C[:CF]=cop2lar(convert(ChainOp,rcopCF))
	end

	if debug_mode
		VIEWCOMPLEX(ret, show=["CV"], explode=[1.2,1.2,1.2], title="arrange3d / ALL atoms")
	end

	ret=SIMPLIFY(ret)
	# @show(ret)
	return ret
end
export ARRANGE3D

# ////////////////////////////////////////////////////////////////
function SPLIT(lar::Lar; debug_mode=false)::Tuple{Lar,Lar}
	if LAR_ARRANGE3D_USE_EXPERIMENTAL
		return arrange3d_v2_split(lar)
	else

		# scrgiorgio: I do not think this is correct, because it could be there is an outer cell with the exact 
		#             same bounding box of an inner cell
		atoms=ATOMS(lar, debug_mode=debug_mode)
		diags =[LinearAlgebra.norm(b[2] - b[1]) for b in [lar_bounding_box(atom, only_used_vertices=true) for atom in atoms]]
		max_diag = maximum(diags)
		outers=lar_copy(lar); outers.C[:CF]=[lar.C[:CF][A] for (A,atom) in enumerate(atoms) if diags[A] == max_diag]
		inners=lar_copy(lar); inners.C[:CF]=[lar.C[:CF][A] for (A,atom) in enumerate(atoms) if diags[A] <  max_diag]
		return outers,inners
	end
end
export SPLIT

# ////////////////////////////////////////////////////////////////
function OUTERS(lar::Lar)::Lar
  return SPLIT(lar)[1]
end
export OUTERS

# ////////////////////////////////////////////////////////////////
function INNERS(lar::Lar)::Lar
  return SPLIT(lar)[2]
end
export INNERS






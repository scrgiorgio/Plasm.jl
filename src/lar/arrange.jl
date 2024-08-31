# //////////////////////////////////////////////////////////////////////////////
function get_external_cycle(V_row::Points, EV::ChainOp, FE::ChainOp)
	FV = abs.(FE) * EV
	vs = sparsevec(mapslices(sum, abs.(EV), dims=1)').nzind
	minv_x1 = maxv_x1 = minv_x2 = maxv_x2 = pop!(vs)
	for i in vs
		if V_row[i, 1] > V_row[maxv_x1, 1]
			maxv_x1 = i
		elseif V_row[i, 1] < V_row[minv_x1, 1]
			minv_x1 = i
		end
		if V_row[i, 2] > V_row[maxv_x2, 2]
			maxv_x2 = i
		elseif V_row[i, 2] < V_row[minv_x2, 2]
			minv_x2 = i
		end
	end
	cells = intersect(
		FV[:, minv_x1].nzind,
		FV[:, maxv_x1].nzind,
		FV[:, minv_x2].nzind,
		FV[:, maxv_x2].nzind
	)
	if length(cells) == 1
		return cells[1]
	else
		for c in cells
			if face_area(V_row, EV, FE[c, :]) < 0
				return c
			end
		end
	end
end

# //////////////////////////////////////////////////////////////////////////////
function minimal_cycles(angles_fn::Function, verbose=true)
	# External interface of TGW algorithm in 2D
	function _minimal_cycles(V::Points,
		ld_bounds::ChainOp)  # , EV)

		lld_cellsnum, ld_cellsnum = size(ld_bounds)
		count_marks = zeros(Int64, ld_cellsnum)
		dir_marks = zeros(Int64, ld_cellsnum)
		d_bounds = spzeros(Int64, ld_cellsnum, 0)

		angles = Array{Array{Int64,1},1}(undef, lld_cellsnum)

		function get_seed_cell()
			s = -1
			for i in 1:ld_cellsnum
				if count_marks[i] == 0
					return i
				elseif count_marks[i] == 1 && s < 0
					s = i
				end
			end
			return s
		end

		for lld in 1:lld_cellsnum
			as = []
			for ld in ld_bounds[lld, :].nzind
				push!(as, (ld, angles_fn(lld, ld)))
			end
			sort!(as, lt=(a, b) -> a[2] < b[2])
			as = map(a -> a[1], as)
			angles[lld] = as
		end

		function nextprev(lld::Int64, ld::Int64, norp)
			as = angles[lld]
			#ne = findfirst(as, ld)  (findfirst(isequal(v), A), 0)[1]
			ne = (findfirst(isequal(ld), as), 0)[1]
			while true
				ne += norp # next or previous
				if ne > length(as)
					ne = 1
				elseif ne < 1
					ne = length(as)
				end

				if count_marks[as[ne]] < 2
					break
				end
			end
			for k = 1:length(count_marks)
				if count_marks[k] > 2
					error("TGW is looping")
				end
			end

			as[ne]
		end

		while (sigma = get_seed_cell()) > 0
			if verbose
				print(Int(floor(50 * sum(count_marks) / ld_cellsnum)), "%\r")
			end

			c_ld = spzeros(Int64, ld_cellsnum)
			if count_marks[sigma] == 0
				c_ld[sigma] = 1
			else
				c_ld[sigma] = -dir_marks[sigma]
			end
			c_lld = ld_bounds * c_ld
			while c_lld.nzind != []
				corolla = spzeros(Int64, ld_cellsnum)
				#corolla = zeros(Int64, ld_cellsnum)

				for tau in c_lld.nzind # when looping, loops here !!
					b_ld = ld_bounds[tau, :]
					pivot = intersect(c_ld.nzind, b_ld.nzind)[1]
					adj = nextprev(tau, pivot, sign(-c_lld[tau]))
					corolla[adj] = c_ld[pivot]
					if b_ld[adj] == b_ld[pivot]
						corolla[adj] *= -1
					end
				end

				c_ld += corolla
				c_lld = ld_bounds * c_ld
			end
			map(s -> count_marks[s] += 1, c_ld.nzind)
			for k = 1:length(count_marks)
				if count_marks[k] > 2
					error("TGW is looping")
				end
			end
			map(s -> dir_marks[s] = c_ld[s], c_ld.nzind)
			d_bounds = [d_bounds c_ld]

		end
		return d_bounds

	end

	return _minimal_cycles
end

# //////////////////////////////////////////////////////////////////////////////
function minimal_2cycles(V::Points, EV::ChainOp)

	function edge_angle(v::Int, e::Int)
		edge = EV[e, :]
		v2 = setdiff(edge.nzind, [v])[1]
		x, y = V[v2, :] - V[v, :]
		return atan(y, x)
	end

	for i in 1:EV.m
		j = min(EV[i, :].nzind...)
		EV[i, j] = -1
	end
	VE = convert(ChainOp, SparseArrays.transpose(EV))
	EF = minimal_cycles(edge_angle)(V, VE)

	return convert(ChainOp, SparseArrays.transpose(EF))
end

# //////////////////////////////////////////////////////////////////////////////
function prune_containment_graph(n, V, EVs, shells, graph)

	for i in 1:n
		an_edge = shells[i].nzind[1]
		origin_index = EVs[i][an_edge, :].nzind[1]
		origin = V[origin_index, :]

		for j in 1:n
			if i != j
				if graph[i, j] == 1
					shell_edge_indexes = shells[j].nzind
					ev = EVs[j][shell_edge_indexes, :]

					if !point_in_face(origin, V, ev)
						graph[i, j] = 0
					end
				end
			end
		end

	end
	return graph
end

# //////////////////////////////////////////////////////////////////////////////
function transitive_reduction!(graph)
	n = size(graph, 1)
	for j in 1:n
		for i in 1:n
			if graph[i, j] > 0
				for k in 1:n
					if graph[j, k] > 0
						graph[i, k] = 0
					end
				end
			end
		end
	end
end

# //////////////////////////////////////////////////////////////////////////////
function pre_containment_test(bboxes)
	n = length(bboxes)
	containment_graph = spzeros(Int8, n, n)

	for i in 1:n
		for j in 1:n
			if i != j && bbox_contains(bboxes[j], bboxes[i])
				containment_graph[i, j] = 1
			end
		end
	end

	return containment_graph
end

# //////////////////////////////////////////////////////////////////////////////
function cell_merging(n, containment_graph, V_row, EVs, boundaries, shells, shell_bboxes)

	V=BYCOL(V_row)

	function bboxes(V_row::Points, indexes::ChainOp)
		boxes = Array{Tuple{Any,Any}}(undef, indexes.n)
		for i in 1:indexes.n
			v_inds = indexes[:, i].nzind
			boxes[i] = boundingbox(V_row[v_inds, :],dims=1)
		end
		boxes
	end
	# initiolization
	sums = Array{Tuple{Int,Int,Int}}(undef, 0)
	# assembling child components with father components
	for father in 1:n
		if sum(containment_graph[:, father]) > 0
			father_bboxes = bboxes(V_row, abs.(EVs[father]') * abs.(boundaries[father]'))
			for child in 1:n
				if containment_graph[child, father] > 0
					child_bbox = shell_bboxes[child]
					for b in 1:length(father_bboxes)
						if bbox_contains(father_bboxes[b], child_bbox)
							push!(sums, (father, b, child))
							break
						end
					end
				end
			end
		end
	end
	# offset assembly initialization
	EV = vcat(EVs...)
	edgenum = size(EV, 1)
	facenum = sum(map(x -> size(x, 1), boundaries))
	FE = spzeros(Int8, facenum, edgenum)
	shells2 = spzeros(Int8, length(shells), edgenum)
	r_offsets = [1]
	c_offset = 1
	# submatrices construction
	for i in 1:n
		min_row = r_offsets[end]
		max_row = r_offsets[end] + size(boundaries[i], 1) - 1
		min_col = c_offset
		max_col = c_offset + size(boundaries[i], 2) - 1
		FE[min_row:max_row, min_col:max_col] = boundaries[i]
		shells2[i, min_col:max_col] = shells[i]
		push!(r_offsets, max_row + 1)
		c_offset = max_col + 1
	end
	# offsetting assembly of component submatrices
	for (f, r, c) in sums
		FE[r_offsets[f]+r-1, :] += shells2[c, :]
	end

	return EV, FE
end

# //////////////////////////////////////////////////////////////////////////////

""" Test of point containment in a polygon face """
function point_in_face(point, V_row::Points, copEV::ChainOp)

	function pointInPolygonClassification(V_row, EV)

		""" Accumulator of partial increments when halfline crosses vertices """
		function crossingTest(new, old, status, count)
			if status == 0
				status = new
				return status, (count + 0.5)
			else
				if status == old
					return 0, (count + 0.5)
				else
					return 0, (count - 0.5)
				end
			end
		end

		""" Set tile code of boxed point w.r.t nine tiles of 2D plane  """
		function setTile(box)
			tiles = [[9, 1, 5], [8, 0, 4], [10, 2, 6]]
			b1, b2, b3, b4 = box
			""" code point position w.r.t query box using Bitwise OR """
			function tileCode(point)
				x, y = point
				code = 0
				if y > b1
					code = code | 1
				end
				if y < b2
					code = code | 2
				end
				if x > b3
					code = code | 4
				end
				if x < b4
					code = code | 8
				end
				return code
			end
			return tileCode
		end

		""" partial function; compute point classification w.r.t polygon edges """
		function pointInPolygonClassification0(pnt)
			x, y = pnt
			xmin, xmax, ymin, ymax = x, x, y, y
			tilecode = setTile([ymax, ymin, xmax, xmin])
			count, status = 0, 0

			for k in 1:EV.m # loop on polygon edges
				edge = EV[k, :]
				p1, p2 = V_row[edge.nzind[1], :], V_row[edge.nzind[2], :]
				(x1, y1), (x2, y2) = p1, p2
				c1, c2 = tilecode(p1), tilecode(p2)
				c_edge, c_un, c_int = xor(c1, c2), c1 | c2, c1 & c2

				if (c_edge == 0) & (c_un == 0)
					return "p_on"
				elseif (c_edge == 12) & (c_un == c_edge)
					return "p_on"
				elseif c_edge == 3
					if c_int == 0
						return "p_on"
					elseif c_int == 4
						count += 1
					end
				elseif c_edge == 15
					x_int = ((y - y2) * (x1 - x2) / (y1 - y2)) + x2
					if x_int > x
						count += 1
					elseif x_int == x
						return "p_on"
					end
				elseif (c_edge == 13) & ((c1 == 4) | (c2 == 4))
					status, count = crossingTest(1, 2, status, count)
				elseif (c_edge == 14) & ((c1 == 4) | (c2 == 4))
					status, count = crossingTest(2, 1, status, count)
				elseif c_edge == 7
					count += 1
				elseif c_edge == 11
					count = count
				elseif c_edge == 1
					if c_int == 0
						return "p_on"
					elseif c_int == 4
						status, count = crossingTest(1, 2, status, count)
					end
				elseif c_edge == 2
					if c_int == 0
						return "p_on"
					elseif c_int == 4
						status, count = crossingTest(2, 1, status, count)
					end
				elseif (c_edge == 4) & (c_un == c_edge)
					return "p_on"
				elseif (c_edge == 8) & (c_un == c_edge)
					return "p_on"
				elseif c_edge == 5
					if (c1 == 0) | (c2 == 0)
						return "p_on"
					else
						status, count = crossingTest(1, 2, status, count)
					end
				elseif c_edge == 6
					if (c1 == 0) | (c2 == 0)
						return "p_on"
					else
						status, count = crossingTest(2, 1, status, count)
					end
				elseif (c_edge == 9) & ((c1 == 0) | (c2 == 0))
					return "p_on"
				elseif (c_edge == 10) & ((c1 == 0) | (c2 == 0))
					return "p_on"
				end
			end
			# final test
			if (round(count) % 2) == 1
				return "p_in"
			else
				return "p_out"
			end
		end
		return pointInPolygonClassification0
	end

	return pointInPolygonClassification(V_row, copEV)(point) == "p_in"
end
export point_in_face

# //////////////////////////////////////////////////////////////////////////////
function delete_edges(edges_to_del, V_row::Points, EV::ChainOp)

	edges_to_keep = setdiff(collect(1:EV.m), edges_to_del)
	EV = EV[edges_to_keep, :]

	all_edges = 1:EV.n
	edges_to_del = Array{Int,1}()
	for i in all_edges
		if length(EV[:, i].nzind) == 0
			push!(edges_to_del, i)
		end
	end

	edges_to_keep = setdiff(all_edges, edges_to_del)
	return V_row[edges_to_keep, :], EV[:, edges_to_keep]
end

# //////////////////////////////////////////////////////////////////////////////
function intersect_edges(V_row::Points, edge1::Chain, edge2::Chain)
	err = 10e-8

	x1, y1, x2, y2 = vcat(map(c -> V_row[c, :], edge1.nzind)...)
	x3, y3, x4, y4 = vcat(map(c -> V_row[c, :], edge2.nzind)...)
	ret = Array{Tuple{Points,Float64},1}()

	v1 = [x2 - x1, y2 - y1]
	v2 = [x4 - x3, y4 - y3]
	v3 = [x3 - x1, y3 - y1]
	ang1 = dot(normalize(v1), normalize(v2))
	ang2 = dot(normalize(v1), normalize(v3))
	parallel = 1 - err < abs(ang1) < 1 + err
	colinear = parallel && (1 - err < abs(ang2) < 1 + err || -err < LinearAlgebra.norm(v3) < err)
	if colinear
		o = [x1 y1]
		v = [x2 y2] - o
		alpha = 1 / dot(v, v')
		ps = [x3 y3; x4 y4]
		for i in 1:2
			a = alpha * dot(v', (reshape(ps[i, :], 1, 2) - o))
			if 0 < a < 1
				push!(ret, (ps[i:i, :], a))
			end
		end
	elseif !parallel
		denom = (v2[2]) * (v1[1]) - (v2[1]) * (v1[2])
		a = ((v2[1]) * (-v3[2]) - (v2[2]) * (-v3[1])) / denom
		b = ((v1[1]) * (-v3[2]) - (v1[2]) * (-v3[1])) / denom

		if -err < a < 1 + err && -err <= b <= 1 + err
			p = [(x1 + a * (x2 - x1)) (y1 + a * (y2 - y1))]
			push!(ret, (p, a))
		end
	end
	return ret
end

# //////////////////////////////////////////////////////////////////////////////
"""
skel_merge(V1::Points, EV1::ChainOp, V2::Points, EV2::ChainOp)

Merge two **1-skeletons**
"""
function skel_merge(V1::Points, EV1::ChainOp, V2::Points, EV2::ChainOp)
	V = [V1; V2]
	EV = blockdiag(EV1, EV2)
	return V, EV
end

"""
skel_merge(V1::Points, EV1::ChainOp, FE1::ChainOp, V2::Points, EV2::ChainOp, FE2::ChainOp)

Merge two **2-skeletons**
"""
function skel_merge(V1::Points, EV1::ChainOp, FE1::ChainOp,V2::Points, EV2::ChainOp, FE2::ChainOp)
	FE = blockdiag(FE1, FE2)
	V, EV = skel_merge(V1, EV1, V2, EV2)
	return V, EV, FE
end

# //////////////////////////////////////////////////////////////////////////////
function buildFV(copEV::ChainOp, face::Chain)
	startv = -1
	nextv = 0
	edge = 0

	vs = Array{Int,1}()

	while startv != nextv
		if startv < 0
			edge = face.nzind[1]
			startv = copEV[edge, :].nzind[face[edge] < 0 ? 2 : 1]
			push!(vs, startv)
		else
			edge = setdiff(intersect(face.nzind, copEV[:, nextv].nzind), edge)[1]
		end
		nextv = copEV[edge, :].nzind[face[edge] < 0 ? 1 : 2]
		push!(vs, nextv)

	end

	return vs[1:end-1]
end
export buildFV

function buildFV(EV::Cells, face::Chain)
	return buildFV(build_copEV(EV), face)
end

# //////////////////////////////////////////////////////////////////////////////
function face_area(V_row::Points, EV::ChainOp, face::Chain)
	function triangle_area(triangle_points::Points)
		ret = ones(3, 3)
		ret[:, 1:2] = triangle_points
		return 0.5 * det(ret)
	end

	area = 0

	fv = buildFV(EV, face)

	verts_num = length(fv)
	v1 = fv[1]

	for i in 2:(verts_num-1)

		v2 = fv[i]
		v3 = fv[i+1]

		area += triangle_area(V_row[[v1, v2, v3], :])
	end

	return area
end

function face_area(V_row::Points, EV::Cells, face::Chain)
	return face_area(V_row, build_copEV(EV), face)
end

# //////////////////////////////////////////////////////////////////////////////
function biconnected_components(EV::ChainOp)

	ps = Array{Tuple{Int,Int,Int},1}()
	es = Array{Tuple{Int,Int},1}()
	todel = Array{Int,1}()
	visited = Array{Int,1}()
	bicon_comps = Array{Array{Int,1},1}()
	hivtx = 1

	function an_edge(point) # TODO: fix bug
		# error? : BoundsError: attempt to access 0×0 SparseMatrix ...
		edges = setdiff(EV[:, point].nzind, todel)
		if length(edges) == 0
			edges = [false]
		end
		edges[1]
	end

	function get_head(edge, tail)
		setdiff(EV[edge, :].nzind, [tail])[1]
	end

	function v_to_vi(v)
		i = findfirst(t -> t[1] == v, ps)
		# seems findfirst changed from 0 to Nothing
		if typeof(i) == Nothing
			return false
		elseif i == 0
			return false
		else
			return ps[i][2]
		end
	end

	push!(ps, (1, 1, 1))
	push!(visited, 1)
	exit = false
	while !exit
		edge = an_edge(ps[end][1])
		if edge != false
			tail = ps[end][2]
			head = get_head(edge, ps[end][1])
			hi = v_to_vi(head)
			if hi == false
				hivtx += 1
				push!(ps, (head, hivtx, ps[end][2]))
				push!(visited, head)
			else
				if hi < ps[end][3]
					ps[end] = (ps[end][1], ps[end][2], hi)
				end
			end
			push!(es, (edge, tail))
			push!(todel, edge)
		else
			if length(ps) == 1
				found = false
				pop!(ps)
				for i in 1:size(EV, 2)
					if !(i in visited)
						hivtx = 1
						push!(ps, (i, hivtx, 1))
						push!(visited, i)
						found = true
						break
					end
				end
				if !found
					exit = true
				end

			else
				if ps[end][3] == ps[end-1][2]
					edges = Array{Int,1}()
					while true
						edge, tail = pop!(es)
						push!(edges, edge)
						if tail == ps[end][3]
							if length(edges) > 1
								push!(bicon_comps, edges)
							end
							break
						end
					end

				else
					if ps[end-1][3] > ps[end][3]
						ps[end-1] = (ps[end-1][1], ps[end-1][2], ps[end][3])
					end
				end
				pop!(ps)
			end
		end
	end
	bicon_comps = sort(bicon_comps, lt=(x, y) -> length(x) > length(y))
	return bicon_comps
end

# //////////////////////////////////////////////////////////////////////////////
function cleandecomposition(V_row, copEV, sigma, edge_map)
	# Deletes edges outside sigma area
	todel = []
	new_edges = []
	map(i -> new_edges = union(new_edges, edge_map[i]), sigma.nzind)
	ev = copEV[new_edges, :]
	for e in 1:copEV.m
		if !(e in new_edges)

			vidxs = copEV[e, :].nzind
			v1, v2 = map(i -> V_row[vidxs[i], :], [1, 2])
			centroid = 0.5 * (v1 + v2)

			if !point_in_face(centroid, V_row, ev)
				push!(todel, e)
			end
		end
	end

	for i in reverse(todel)
		for row in edge_map
			filter!(x -> x != i, row)
			for j in 1:length(row)
				if row[j] > i
					row[j] -= 1
				end
			end
		end
	end

	V_row, copEV = delete_edges(todel, V_row, copEV)
	return V_row, copEV
end

# //////////////////////////////////////////////////////////////////////////////
using NearestNeighbors

function merge_vertices!(V::Points, EV::ChainOp, edge_map, err=1e-4)
	vertsnum = size(V, 1)
	edgenum = size(EV, 1)
	newverts = zeros(Int, vertsnum)

	# NearestNeighbors.KDTree constructor needs an explicit array of Float64
	V = Array{Float64,2}(V)
	kdtree = NearestNeighbors.KDTree(BYROW(V))

	# merge congruent vertices
	todelete = []
	i = 1
	for vi in 1:vertsnum
		if !(vi in todelete)
			nearvs = NearestNeighbors.inrange(kdtree, V[vi, :], err)
			newverts[nearvs] .= i
			nearvs = setdiff(nearvs, vi)
			todelete = union(todelete, nearvs)
			i = i + 1
		end
	end
	nV = V[setdiff(collect(1:vertsnum), todelete), :]

	# merge congruent edges
	edges = Array{Tuple{Int,Int},1}(undef, edgenum)
	oedges = Array{Tuple{Int,Int},1}(undef, edgenum)
	for ei in 1:edgenum
		v1, v2 = EV[ei, :].nzind
		edges[ei] = Tuple{Int,Int}(sort([newverts[v1], newverts[v2]]))
		oedges[ei] = Tuple{Int,Int}(sort([v1, v2]))
	end
	nedges = union(edges)
	nedges = filter(t -> t[1] != t[2], nedges)
	nedgenum = length(nedges)
	nEV = spzeros(Int8, nedgenum, size(nV, 1))
	# maps pairs of vertex indices to edge index
	etuple2idx = Dict{Tuple{Int,Int},Int}()
	# builds `edge_map`
	for ei in 1:nedgenum
		nEV[ei, collect(nedges[ei])] .= 1
		etuple2idx[nedges[ei]] = ei
	end
	for i in 1:length(edge_map)
		row = edge_map[i]
		row = map(x -> edges[x], row)
		row = filter(t -> t[1] != t[2], row)
		row = map(x -> etuple2idx[x], row)
		edge_map[i] = row
	end
	# return new vertices and new edges
	return Points(nV), nEV
end

# //////////////////////////////////////////////////////////////////////////////
function frag_edge(V_row, EV::ChainOp, edge_idx::Int, bigPI)
	alphas = Dict{Float64,Int}()
	edge = EV[edge_idx, :]
	verts = V_row[edge.nzind, :]
	for i in bigPI[edge_idx]
		if i != edge_idx
			intersection = intersect_edges(V_row, edge, EV[i, :])
			for (point, alpha) in intersection
				verts = [verts; point]
				alphas[alpha] = size(verts, 1)
			end
		end
	end
	alphas[0.0], alphas[1.0] = [1, 2]
	alphas_keys = sort(collect(keys(alphas)))
	edge_num = length(alphas_keys) - 1
	verts_num = size(verts, 1)
	ev = SparseArrays.spzeros(Int8, edge_num, verts_num)
	for i in 1:edge_num
		ev[i, alphas[alphas_keys[i]]] = 1
		ev[i, alphas[alphas_keys[i+1]]] = 1
	end
	return verts, ev
end

# //////////////////////////////////////////////////////////////////////////////
function permutationOrbits(perm::OrderedDict)
	out = Array{Int64,1}[]
	while perm ≠ Dict()
		x = collect(keys(perm))[1]
		orbit = Int64[]
		while x in keys(perm)
			append!(orbit, perm[x])
			y, x = x, perm[x]
			delete!(perm, y)
		end
		append!(out, [push!(orbit, orbit[1])])
	end
	return out
end

# //////////////////////////////////////////////////////////////////////////////
function faces2polygons(copEV, copFE)
	polygons = Array{Array{Int64,1},1}[]
	cycles = Array{Array{Array{Int64,1},1},1}[]
	for f = 1:size(copFE, 1)
		edges, signs = findnz(copFE[f, :])
		permutationMap = OrderedDict([s > 0 ? findnz(copEV[e, :])[1] : reverse(findnz(copEV[e, :])[1])
																	for (e, s) in zip(edges, signs)])
		orbits = permutationOrbits(permutationMap)
		edgecycles = [[[orbit[k], orbit[k+1]] for k = 1:length(orbit)-1] for orbit in orbits]
		push!(polygons, [orbit[1:end-1] for orbit in orbits])
		push!(cycles, edgecycles)
	end
	return polygons, cycles
end


# //////////////////////////////////////////////////////////////////////////////
function planar_arrangement(V::Points, copEV::ChainOp, sigma::Chain=spzeros(Int8, 0))

	V_row = BYROW(V)

	# planar_arrangement_1

	# data structures initialization
	edgenum = size(copEV, 1)
	edge_map = Array{Array{Int,1},1}(undef, edgenum)
	rV = Points(zeros(0, 2))
	rEV = SparseArrays.spzeros(Int8, 0, 0)
	finalcells_num = 0

	# space index computation
	bigPI = spaceindex(V_row, cop2lar(copEV))

	# sequential (iterative) processing of edge fragmentation
	for i in 1:edgenum
		v, ev = frag_edge(V_row, copEV, i, bigPI)
		newedges_nums = map(x -> x + finalcells_num, collect(1:size(ev, 1)))
		edge_map[i] = newedges_nums
		finalcells_num += size(ev, 1)
		rV = convert(Points, rV)
		rV, rEV = skel_merge(rV, rEV, v, ev)
	end
	#  end
	# merging of close vertices and edges (2D congruence)
	V_row, copEV = rV, rEV
	V_row, copEV = merge_vertices!(V_row, copEV, edge_map)

	# cleandecomposition
	if sigma.n > 0
		V_row, copEV = cleandecomposition(V_row, copEV, sigma, edge_map)
	end

	bicon_comps = biconnected_components(copEV)

	if isempty(bicon_comps)
		return (nothing, nothing, nothing)
	end

	function planar_arrangement_2(V_row, copEV, bicon_comps, edge_map, sigma::Chain=spzeros(Int8, 0))

		edges = sort(union(bicon_comps...))
		todel = sort(setdiff(collect(1:size(copEV, 1)), edges))

		for i in reverse(todel)
			for row in edge_map

				filter!(x -> x != i, row)

				for j in 1:length(row)
					if row[j] > i
						row[j] -= 1
					end
				end
			end
		end


		bicon_comps = biconnected_components(copEV)

		function componentgraph(V_row, copEV, bicon_comps)
			# arrangement of isolated components
			n = size(bicon_comps, 1)
			shells = Array{Chain,1}(undef, n)
			boundaries = Array{ChainOp,1}(undef, n)
			EVs = Array{ChainOp,1}(undef, n)
			# for each component
			for p in 1:n
				ev = copEV[sort(bicon_comps[p]), :]
				# computation of 2-cells
				fe = minimal_2cycles(V_row, ev)
				# exterior cycle
				shell_num = get_external_cycle(V_row, ev, fe)
				# decompose each fe (co-boundary local to component)
				EVs[p] = ev
				tokeep = setdiff(1:fe.m, shell_num)
				boundaries[p] = fe[tokeep, :]
				shells[p] = fe[shell_num, :]
			end
			# computation of bounding boxes of isolated components
			shell_bboxes = []
			for i in 1:n
				vs_indexes = (abs.(EVs[i]') * abs.(shells[i])).nzind
				push!(shell_bboxes, boundingbox(V_row[vs_indexes, :],dims=1))
			end
			# computation and reduction of containment graph
			containment_graph = pre_containment_test(shell_bboxes)
			containment_graph = prune_containment_graph(n, V_row, EVs, shells, containment_graph)
			transitive_reduction!(containment_graph)
			return n, containment_graph, V_row, EVs, boundaries, shells, shell_bboxes
		end

		n, containment_graph, V_row, EVs, boundaries, shells, shell_bboxes = componentgraph(V_row, copEV, bicon_comps)

		copEV, FE = cell_merging(n, containment_graph, V_row, EVs, boundaries, shells, shell_bboxes)

		return BYCOL(V_row), copEV, FE
	end


	#Planar_arrangement_2
	return planar_arrangement_2(V_row, copEV, bicon_comps, edge_map, sigma)
end


# //////////////////////////////////////////////////////////////////////////////
function arrange2D(V, EV)

	copEV = coboundary_0(EV::Cells)
	cop_EW = convert(ChainOp, copEV)
	V, copEV, copFE = planar_arrangement(V, cop_EW::ChainOp)
	EVs = FV2EVs(copEV, copFE) # polygonal face fragments

	V_row = BYROW(V)

	# triangulate
	triangles_per_face = Array{Any,1}(undef, copFE.m)
	if size(V_row, 2) == 2
		V_row = [V_row zeros(size(V_row, 1), 1)]
	end

	polygons, edgecycles = faces2polygons(copEV, copFE) #new

	for f in 1:copFE.m
		edges_idxs = copFE[f, :].nzind
		edge_num = length(edges_idxs)
		edges = Array{Int,1}[] #zeros(Int, edge_num, 2)

		# fv = buildFV(copEV, copFE[f, :])
		fv = union(polygons[f]...)
		vs = V_row[fv, :]
		edges = union(edgecycles[f]...)
		edges = convert(Array{Int,2}, hcat(edges...)')

		v = convert(Points, vs'[1:2, :])
		vmap = Dict(zip(fv, 1:length(fv))) # vertex map
		mapv = Dict(zip(1:length(fv), fv)) # inverse vertex map
		ev = [[vmap[e] for e in edges[k, :]] for k = 1:size(edges, 1)]
		triangles_per_face[f] = [[mapv[jt] for jt in triangle] for triangle in TRIANGULATE2D(v, ev)]

		tV = V_row[:, 1:2]

		area = face_area(tV, copEV, copFE[f, :])
		if area < 0
			for i in 1:length(triangles_per_face[f])
				triangles_per_face[f][i] = triangles_per_face[f][i][end:-1:1]
			end
		end
	end

	FVs = convert(Array{Cells}, triangles_per_face)
	return V, FVs, EVs, copEV, copFE
end
export arrange2D

# //////////////////////////////////////////////////////////////////////////////
function u_coboundary_1(copFV::ChainOp, copEV::ChainOp, convex=true::Bool)::ChainOp
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
	return copFE
end

# //////////////////////////////////////////////////////////////////////////////
function u_coboundary_1(FV::Cells, EV::Cells, convex=true::Bool)::ChainOp
	copFV = lar2cop(FV)
	copEV = lar2cop(EV)
	out = u_coboundary_1(copFV::ChainOp, copEV::ChainOp, convex::Bool)
	return out
end

# //////////////////////////////////////////////////////////////////////////////
export coboundary_1
function coboundary_1(V::Points, copFV::ChainOp, copEV::ChainOp, convex=true::Bool, exterior=false::Bool)::ChainOp

	copFE = u_coboundary_1(copFV::ChainOp, copEV::ChainOp, convex)
	EV = [findnz(copEV[k, :])[1] for k = 1:size(copEV, 1)]
	copEV = sparse(coboundary_0(EV::Cells))
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
		EV = convert(ChainOp, SparseArrays.transpose(boundary_1(EV)))

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


function coboundary_1(V::Points, FV::Cells, EV::Cells; convex=true::Bool, exterior=false::Bool)::ChainOp
	# generate unsigned operator's sparse matrix
	copFV = lar2cop(FV)
	copEV = lar2cop(EV)
	# greedy generation of incidence signs
	return coboundary_1(V, copFV::ChainOp, copEV::ChainOp, convex, exterior)
end

# //////////////////////////////////////////////////////////////////////////////
""" Main function of arrangement pipeline """
function space_arrangement(V::Points, EV::ChainOp, FE::ChainOp)

	function frag_face(V, EV::ChainOp, FE::ChainOp, sp_idx, sigma)
		vs_num = size(V, 1)
		# 2D transformation of `sigma` face
		sigmavs = (abs.(FE[sigma:sigma, :])*abs.(EV))[1, :].nzind # sigma vertex indices
		sV = V[sigmavs, :]
		sEV = EV[FE[sigma, :].nzind, sigmavs]
		M = submanifold_mapping(sV)
		tV = ([V ones(vs_num)]*M)[:, 1:3]  # folle convertire *tutti* i vertici
		sV = tV[sigmavs, :]
		# `sigma` face intersection with faces in `sp_idx[sigma]`, i.e., in `bigpi`
		for i in sp_idx[sigma]
			tmpV, tmpEV = face_int(tV, EV, FE[i, :]) # va in loop qui dentro
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

	# historically arrangement works internally by using by-row vertices
	V_row = BYROW(V)
	fs_num = size(FE, 1)
	# strange but necessary cycle of computations to get FV::Cells algebraically
	FV = (abs.(FE) * abs.(EV)) .÷ 2
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
		nV, nEV, nFE = frag_face(V_row, EV, FE, sp_idx, sigma)
		depot_V[sigma] = nV
		depot_EV[sigma] = nEV
		depot_FE[sigma] = nFE
	end
	rV = vcat(depot_V...)
	rEV = SparseArrays.blockdiag(depot_EV...)
	rFE = SparseArrays.blockdiag(depot_FE...)
	rV, rcopEV, rcopFE = Plasm.merge_vertices(rV, rEV, rFE)
	rcopCF = build_copFC(rV, rcopEV, rcopFE)

	# historically arrangement works internally by using by-row vertices
	return BYCOL(rV), rcopEV, rcopFE, rcopCF
end
export space_arrangement

# //////////////////////////////////////////////////////////////////////////////
function spaceindex(V_row::Points, CV)::Cells
	V = BYCOL(V_row)
	pdim = size(V, 1)
	cellpoints = [V[:, CV[k]] for k = 1:LEN(CV)]
	bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints]
	boxdicts = [bbox_coord_intervals(k, bboxes) for k = 1:pdim]
	trees = [IntervalTrees.IntervalMap{Float64,Array}() for k = 1:pdim]
	covers = []
	for k = 1:pdim
		for (key, boxset) in boxdicts[k]
			trees[k][tuple(key...)] = boxset
		end
		push!(covers, bbox_covering(bboxes, k, trees[k]))
	end
	covers = [reduce(intersect, [covers[h][k] for h = 1:pdim]) for k = 1:length(CV)]
end


#//////////////////////////////////////////////////////////////////////////////
""" Compute the map from vs (at least three) to z=0 """
# REMARK: will not work when vs[:,1:3] are aligned !!!!  TODO: fix 
function submanifold_mapping(vs)
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
function face_int(V::Points, EV::ChainOp, face::Chain)
	vs = buildFV(EV, face)     # EV::ChainOp, face::Chain
	retV = Points(undef, 0, 3)
	err = 10e-8
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
""" Task to iteratively add new local components to the global 2-skeleton """
# Remark: sensible to `err`; works w `err=1e-6`, "ERROR: no pivot" with `err=1e-7`
function merge_vertices(V::Points, EV::ChainOp, FE::ChainOp, err=1e-6)
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
""" Component of TGW in 3D (Pao);  return an ordered `fan` of 2-cells """
function ord(hinge::Int, bd1, V::Points, FV::Cells, EV::Cells, FE::Cells)

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
			if interior_to_f([v1, v2, v3], f, V, FV, EV, FE)
				push!(triangles, [v1, v2, v3])
				break
			else
				index = findnext(v -> (area(v1, v2, v) ≠ 0), FV[f], index + 1)
				v3 = FV[f][index]
			end
		end
	end
	order = ordering(triangles, V)
	return [cells[index] for index in order]
end


# //////////////////////////////////////////////////////////////////////////////
""" TGW algorithm implementation (pao) """
function build_copFC(V_row, rcopEV, rcopFE)

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
				fan = ord(abs(τ), bd1, V, FV, EV, FE) # ord(pivot,bd1)
				if τ > 0
					adj = mynext(fan, pivot, marks)
				elseif τ < 0
					adj = myprev(fan, pivot, marks)
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
	return copCF
end


# //////////////////////////////////////////////////////////////////////////////
""" Triangle containment test of checkpoint into `f`; used in TGW algorithm """
function interior_to_f(triangle, f, V, FV, EV, FE)::Bool
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
""" Correct ordering of triangles (hence faces) about each boundary edge """
function ordering(triangles, V)
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
""" Utility function for TGW in 3D """
function mynext(cycle, pivot, marks)
	len = length(cycle)
	ind = findfirst(x -> x == pivot, cycle)[1]
	nextIndex = ind == len ? 1 : ind + 1
	#if marks[nextIndex]==2
	#  nextIndex = ind==1 ? len : ind-1
	#end
	return cycle[nextIndex][1]
end

# //////////////////////////////////////////////////////////////////////////////
""" Utility function for TGW in 3D """
function myprev(cycle, pivot, marks)
	len = length(cycle)
	ind = findfirst(x -> x == pivot, cycle)[1] #  ind is address of pivot in cycle
	nextIndex = ind == 1 ? len : ind - 1
	#if marks[nextIndex]==2
	#  nextIndex = ind==len ? 1 : ind+1
	#end
	return cycle[nextIndex][1]
end

# //////////////////////////////////////////////////////////////////////////////
""" Coherently orient the edges of f face """
function vcycle(copEV::ChainOp, copFE::ChainOp, f::Int)
	edges, signs = findnz(copFE[f, :])
	vpairs = [s > 0 ? findnz(copEV[e, :])[1] :
						reverse(findnz(copEV[e, :])[1])
						for (e, s) in zip(edges, signs)]
	a = [pair for pair in vpairs if length(pair) == 2]
	function mycat(a::Cells)
		out = []
		for cell in a
			append!(out, cell)
		end
		return out
	end
	vs = collect(Set(mycat(a)))
	vdict = Dict(zip(vs, 1:length(vs)))
	edges = [[vdict[pair[1]], vdict[pair[2]]] for pair in vpairs if length(pair) == 2]
	return vs, edges
end
export vcycle




using IntervalTrees

# ///////////////////////////////////////////////////////////////////
"""
lar_bbox_create(vertices::Points)

The axis aligned *bounding box* of the provided Matrix of n-dim `vertices`.
The box is returned as the pair of `Points` of two opposite corners.

NOTE: assuming LAR by-col representation, if it's by-row using dims=1
"""
function lar_bbox_create(vertices::Points;dims::Int=2)
	minimum = mapslices(x -> min(x...), vertices, dims=dims)
	maximum = mapslices(x -> max(x...), vertices, dims=dims)
	return minimum, maximum
end
export lar_bbox_create

function lar_bbox_contains(container, contained)
	b1_min, b1_max = container
	b2_min, b2_max = contained
	all(map((i, j, k, l) -> i <= j <= k <= l, b1_min, b2_min, b2_max, b1_max))
end
export lar_bbox_contains

""" E' utilizzata da function spaceindex(V_row::Points, CV)::Cells 
NOTE:  CV è generico: si può usare per C, F, E ... calcola covers (lista delle celle i cui box possono coprire quelli  in bboxes)
"""
function lar_bbox_covering(bboxes, index, tree)
	covers = [[] for k = 1:length(bboxes)] #ata da apce ??
	for (i, boundingbox) in enumerate(bboxes)
		extent = bboxes[i][index, :]
		iterator = IntervalTrees.intersect(tree, tuple(extent...))
		for x in iterator
			append!(covers[i], x.value)
		end
	end
	return covers
end
export lar_bbox_covering

""" Make dictionary of 1D boxes for IntervalTrees construction """
function lar_bbox_coord_intervals(coord, bboxes)
	boxdict = OrderedDict{Array{Float64,1},Array{Int64,1}}()
	for (h, box) in enumerate(bboxes)
		key = box[coord, :]
		if haskey(boxdict, key) == false
			boxdict[key] = [h]
		else
			push!(boxdict[key], h)
		end
	end
	return boxdict
end
export lar_bbox_coord_intervals

# //////////////////////////////////////////////////////////////////////////////
function lar_bbox_containment_graph(bboxes)
	n = length(bboxes)
	ret = spzeros(Int8, n, n)
	for i in 1:n
		for j in 1:n
			if i != j && lar_bbox_contains(bboxes[j], bboxes[i])
				ret[i, j] = 1
			end
		end
	end
	return ret
end

# //////////////////////////////////////////////////////////////////////////////
function spaceindex(V_row::Points, CV)::Cells
	V = BYCOL(V_row)
	pdim = size(V, 1)
	cellpoints = [V[:, CV[k]] for k = 1:LEN(CV)]
	bboxes = [hcat(lar_bbox_create(cell)...) for cell in cellpoints]
	boxdicts = [lar_bbox_coord_intervals(k, bboxes) for k = 1:pdim]
	trees = [IntervalTrees.IntervalMap{Float64,Array}() for k = 1:pdim]
	covers = []
	for k = 1:pdim
		for (key, boxset) in boxdicts[k]
			trees[k][tuple(key...)] = boxset
		end
		push!(covers, lar_bbox_covering(bboxes, k, trees[k]))
	end
	covers = [reduce(intersect, [covers[h][k] for h = 1:pdim]) for k = 1:length(CV)]
end


# //////////////////////////////////////////////////////////////////////////////
function intersect_edges(V_row::Points, edge1::Chain, edge2::Chain; err=LAR_DEFAULT_ERR)

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
function merge_vertices_2d!(V::Points, EV::ChainOp, edge_map; err=LAR_DEFAULT_ERR)
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
function clean_decomposition(V_row, copEV, sigma, edge_map)
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
	function _minimal_cycles(V::Points, ld_bounds::ChainOp)  # , EV)

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
function minimal_2d_cycles(V::Points, EV::ChainOp)

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
function cell_merging(n, containment_graph, V_row, EVs, boundaries, shells, shell_bboxes)

	V=BYCOL(V_row)

	function bboxes(V_row::Points, indexes::ChainOp)
		boxes = Array{Tuple{Any,Any}}(undef, indexes.n)
		for i in 1:indexes.n
			v_inds = indexes[:, i].nzind
			boxes[i] = lar_bbox_create(V_row[v_inds, :],dims=1)
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
						if lar_bbox_contains(father_bboxes[b], child_bbox)
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
"""
skel_merge(V1::Points, EV1::ChainOp, V2_row::Points, EV2::ChainOp)

Merge two **1-skeletons**
"""
function skel_merge(V1_row::Points, EV1::ChainOp, V2_row::Points, EV2::ChainOp)
	V_row = [V1_row; V2_row]
	EV = blockdiag(EV1, EV2)
	return V_row, EV
end

"""
skel_merge(V1_row::Points, EV1::ChainOp, FE1::ChainOp, V2::Points, EV2::ChainOp, FE2::ChainOp)

Merge two **2-skeletons**
"""
function skel_merge(V1_row::Points, EV1::ChainOp, FE1::ChainOp,V2_row::Points, EV2::ChainOp, FE2::ChainOp)
	FE = blockdiag(FE1, FE2)
	V_row, EV = skel_merge(V1_row, EV1, V2_row, EV2)
	return V_row, EV, FE
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
	V_row, copEV = merge_vertices_2d!(V_row, copEV, edge_map)

	if sigma.n > 0
		V_row, copEV = clean_decomposition(V_row, copEV, sigma, edge_map)
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
				fe = minimal_2d_cycles(V_row, ev)
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
				push!(shell_bboxes, lar_bbox_create(V_row[vs_indexes, :],dims=1))
			end
			# computation and reduction of containment graph
			containment_graph = lar_bbox_containment_graph(shell_bboxes)
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


# //////////////////////////////////////////////////////////////////////////////
function ARRANGE2D(lar::Lar)
	copEV = convert(ChainOp, cop_coboundary_0(lar.C[:EV]))
	V, copEV, copFE = planar_arrangement(lar.V, copEV)
	ret=Lar(V, Dict( 
		:EV => cop2lar(copEV), 
		:FE => cop2lar(copFE)))
	ret.C[:FV]=compute_FV(ret)
	return SIMPLIFY(ret)

end
export ARRANGE2D



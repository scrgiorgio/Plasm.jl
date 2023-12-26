using LinearAlgebra
using SparseArrays
using DataStructures
using NearestNeighbors
using Triangulate
using IntervalTrees

export arrange2D, LAR, Points, Cells, Cell, Chain, ChainOp, ChainComplex,
bbox

const Points = Matrix
const Cells = Vector{Vector{Int}}
const Cell = SparseVector{Int8, Int}
const Chain = SparseVector{Int8,Int}
const ChainOp = SparseMatrixCSC{Int8,Int}
const ChainComplex = Vector{ChainOp}


# //////////////////////////////////////////////////////////////////////////////
"""
    bbox(vertices::Points)

The axis aligned *bounding box* of the provided Matrix of n-dim `vertices`.

The box is returned as the pair of `Points` of two opposite corners.
"""
function bbox(vertices::Points)
    minimum = mapslices(x->min(x...), vertices, dims=1)
    maximum = mapslices(x->max(x...), vertices, dims=1)
    minimum, maximum
end

# //////////////////////////////////////////////////////////////////////////////

""" 
   get_external_cycle(V::Points, EV::ChainOp, FE::ChainOp)
   
Find external cycles (area < 0) in chain complex (FE already computed).
"""
function get_external_cycle(V::Points, EV::ChainOp, FE::ChainOp)
println("get_external_cycle")
    FV = abs.(FE)*EV
    vs = sparsevec(mapslices(sum, abs.(EV), dims=1)').nzind
    minv_x1 = maxv_x1 = minv_x2 = maxv_x2 = pop!(vs)
    for i in vs
        if V[i, 1] > V[maxv_x1, 1]
            maxv_x1 = i
        elseif V[i, 1] < V[minv_x1, 1]
            minv_x1 = i
        end
        if V[i, 2] > V[maxv_x2, 2]
            maxv_x2 = i
        elseif V[i, 2] < V[minv_x2, 2]
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
            if face_area(V, EV, FE[c, :]) < 0
                return c
            end
        end
    end
end

# //////////////////////////////////////////////////////////////////////////////
"""
   minimal_cycles(V::Points, EV::ChainOp)(V::Points, EV::ChainOp)
   
Interface of TGW algorithm in 2D
"""
function minimal_cycles(angles_fn::Function, verbose=true)
println("minimal_cycles")
    # External interface of TGW algorithm in 2D
    function _minimal_cycles(V::Points,
    	  ld_bounds::ChainOp)  # , EV)
println("_minimal_cycles")

        lld_cellsnum, ld_cellsnum = size(ld_bounds)
        count_marks = zeros(Int64, ld_cellsnum)
        dir_marks = zeros(Int64, ld_cellsnum)
        d_bounds = spzeros(Int64, ld_cellsnum, 0)

        angles = Array{Array{Int64, 1}, 1}(undef, lld_cellsnum)

        function get_seed_cell()
println("get_seed_cell")
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
            sort!(as, lt=(a,b)->a[2]<b[2])
            as = map(a->a[1], as)
            angles[lld] = as
        end

        function nextprev(lld::Int64, ld::Int64, norp)
println("nextprev")
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
            for k=1:length(count_marks)
            	if count_marks[k]>2  error("TGW is looping") end
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
            c_lld = ld_bounds*c_ld
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
                c_lld = ld_bounds*c_ld
            end
            map(s->count_marks[s] += 1, c_ld.nzind)
            for k=1:length(count_marks)
            	if count_marks[k]>2  error("TGW is looping") end
            end
            map(s->dir_marks[s] = c_ld[s], c_ld.nzind)
            d_bounds = [d_bounds c_ld]

        end
        return d_bounds

    end

    return _minimal_cycles
end

# //////////////////////////////////////////////////////////////////////////////

function minimal_2cycles(V::Points, EV::ChainOp)
println("minimal_2cycles")

    function edge_angle(v::Int, e::Int)
        edge = EV[e, :]
        v2 = setdiff(edge.nzind, [v])[1]
        x, y = V[v2, :] - V[v, :]
        return atan(y, x)
    end

    for i in 1:EV.m
        j = min(EV[i,:].nzind...)
        EV[i, j] = -1
    end
    VE = convert(ChainOp, SparseArrays.transpose(EV))
    EF = minimal_cycles(edge_angle)(V, VE)

    return convert(ChainOp, SparseArrays.transpose(EF))
end

# //////////////////////////////////////////////////////////////////////////////

function bbox_contains(container, contained)
println("bbox_contains")
    b1_min, b1_max = container
    b2_min, b2_max = contained
    all(map((i,j,k,l)->i<=j<=k<=l, b1_min, b2_min, b2_max, b1_max))
end

# //////////////////////////////////////////////////////////////////////////////

function prune_containment_graph(n, V, EVs, shells, graph)
println("prune_containment_graph")

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
println("transitive_reduction")
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
println("pre_containment_test")
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

function componentgraph(V, copEV, bicon_comps)
println("componentgraph")
    # arrangement of isolated components
	n = size(bicon_comps, 1)
   	shells = Array{Chain, 1}(undef, n)
	boundaries = Array{ChainOp, 1}(undef, n)
	EVs = Array{ChainOp, 1}(undef, n)
    # for each component
	for p in 1:n
		ev = copEV[sort(bicon_comps[p]), :]
        # computation of 2-cells
		fe = minimal_2cycles(V, ev)
        # exterior cycle
		shell_num = get_external_cycle(V, ev, fe)
        # decompose each fe (co-boundary local to component)
		EVs[p] = ev
		tokeep = setdiff(1:fe.m, shell_num)
		boundaries[p] = fe[tokeep, :]
		shells[p] = fe[shell_num, :]
    end
    # computation of bounding boxes of isolated components
	shell_bboxes = []
	for i in 1:n
    	vs_indexes = (abs.(EVs[i]')*abs.(shells[i])).nzind
   		push!(shell_bboxes, bbox(V[vs_indexes, :]))
	end
    # computation and reduction of containment graph
	containment_graph = pre_containment_test(shell_bboxes)
	containment_graph = prune_containment_graph(n, V, EVs, shells, containment_graph)
	transitive_reduction!(containment_graph)
	return n, containment_graph, V, EVs, boundaries, shells, shell_bboxes
end

# //////////////////////////////////////////////////////////////////////////////

function cell_merging(n,containment_graph,V,EVs,boundaries,shells,shell_bboxes)
println("cell_merging")
    function bboxes(V::Points, indexes::ChainOp)
println("bboxes")
        boxes = Array{Tuple{Any, Any}}(undef, indexes.n)
        for i in 1:indexes.n
            v_inds = indexes[:, i].nzind
            boxes[i] = bbox(V[v_inds, :])
        end
        boxes
    end
    # initiolization
    sums = Array{Tuple{Int, Int, Int}}(undef, 0);
    # assembling child components with father components
    for father in 1:n
        if sum(containment_graph[:, father]) > 0
            father_bboxes = bboxes(V, abs.(EVs[father]')*abs.(boundaries[father]'))
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
    facenum = sum(map(x->size(x,1), boundaries))
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

function boundingbox(vertices::Points)
println("boundingbox")
   minimum = mapslices(x->min(x...), vertices, dims=2)
   maximum = mapslices(x->max(x...), vertices, dims=2)
   return minimum, maximum
end

# //////////////////////////////////////////////////////////////////////////////

function coordintervals(coord,bboxes)
println("boundingbox")
	boxdict = OrderedDict{Array{Float64,1},Array{Int64,1}}()
	for (h,box) in enumerate(bboxes)
		key = box[coord,:]
		if haskey(boxdict,key) == false
			boxdict[key] = [h]
		else
			push!(boxdict[key], h)
		end
	end
	return boxdict
end

# //////////////////////////////////////////////////////////////////////////////

function boxcovering(bboxes, index, tree)
println("boxcovering")
  covers = [[] for k=1:length(bboxes)]
  for (i, boundingbox) in enumerate(bboxes)
	extent = bboxes[i][index,:]
	iterator = IntervalTrees.intersect(tree, tuple(extent...))
	for x in iterator
	  append!(covers[i], x.value)
	end
  end
  return covers
end

# //////////////////////////////////////////////////////////////////////////////

function lar2cop(CV::Cells)::ChainOp
println("boxcovering")
	I = Int[]; J = Int[]; Value = Int8[];
	for k=1:size(CV,1)
		n = length(CV[k])
		append!(I, k * ones(Int, n))
		append!(J, CV[k])
		append!(Value, ones(Int, n))
	end
	return SparseArrays.sparse(I,J,Value)
end

# //////////////////////////////////////////////////////////////////////////////

function point_in_face(point, V::Points, copEV::ChainOp)
println("point_in_face")

    function pointInPolygonClassification(V,EV)
println("pointInPolygonClassification")

        function crossingTest(new, old, status, count)
println("crossingTest")
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

        function setTile(box)
println("setTile")
        tiles = [[9,1,5],[8,0,4],[10,2,6]]
        b1,b2,b3,b4 = box
        function tileCode(point)
println("tileCode")
            x,y = point
            code = 0
            if y>b1 code=code|1 end
            if y<b2 code=code|2 end
            if x>b3 code=code|4 end
            if x<b4 code=code|8 end
            return code
        end
        return tileCode
        end

        function pointInPolygonClassification0(pnt)
println("pointInPolygonClassification0")
            x,y = pnt
            xmin,xmax,ymin,ymax = x,x,y,y
            tilecode = setTile([ymax,ymin,xmax,xmin])
            count,status = 0,0

            for k in 1:EV.m
                edge = EV[k,:]
                p1, p2 = V[edge.nzind[1], :], V[edge.nzind[2], :]
                (x1,y1),(x2,y2) = p1,p2
                c1,c2 = tilecode(p1),tilecode(p2)
                c_edge, c_un, c_int = xor(c1, c2), c1|c2, c1&c2

                if (c_edge == 0) & (c_un == 0) return "p_on"
                elseif (c_edge == 12) & (c_un == c_edge) return "p_on"
                elseif c_edge == 3
                    if c_int == 0 return "p_on"
                    elseif c_int == 4 count += 1 end
                elseif c_edge == 15
                    x_int = ((y-y2)*(x1-x2)/(y1-y2))+x2
                    if x_int > x count += 1
                    elseif x_int == x return "p_on" end
                elseif (c_edge == 13) & ((c1==4) | (c2==4))
                        status, count = crossingTest(1,2,status,count)
                elseif (c_edge == 14) & ((c1==4) | (c2==4))
                        status, count = crossingTest(2,1,status,count)
                elseif c_edge == 7 count += 1
                elseif c_edge == 11 count = count
                elseif c_edge == 1
                    if c_int == 0 return "p_on"
                    elseif c_int == 4
                        status, count = crossingTest(1,2,status,count)
                    end
                elseif c_edge == 2
                    if c_int == 0 return "p_on"
                    elseif c_int == 4
                        status, count = crossingTest(2,1,status,count)
                    end
                elseif (c_edge == 4) & (c_un == c_edge) return "p_on"
                elseif (c_edge == 8) & (c_un == c_edge) return "p_on"
                elseif c_edge == 5
                    if (c1==0) | (c2==0) return "p_on"
                    else
                        status, count = crossingTest(1,2,status,count)
                    end
                elseif c_edge == 6
                    if (c1==0) | (c2==0) return "p_on"
                    else
                        status, count = crossingTest(2,1,status,count)
                    end
                elseif (c_edge == 9) & ((c1==0) | (c2==0)) return "p_on"
                elseif (c_edge == 10) & ((c1==0) | (c2==0)) return "p_on"
                end
            end

            if (round(count)%2)==1
                return "p_in"
            else
                return "p_out"
            end
        end
        return pointInPolygonClassification0
    end

    return pointInPolygonClassification(V, copEV)(point) == "p_in"
end

# //////////////////////////////////////////////////////////////////////////////

function delete_edges(todel, V::Points, EV::ChainOp)
println("delete_edges")
    tokeep = setdiff(collect(1:EV.m), todel)
    EV = EV[tokeep, :]

    vertinds = 1:EV.n
    todel = Array{Int, 1}()
    for i in vertinds
        if length(EV[:, i].nzind) == 0
            push!(todel, i)
        end
    end

    tokeep = setdiff(vertinds, todel)
    EV = EV[:, tokeep]
    V = V[tokeep, :]

    return V, EV
end

# //////////////////////////////////////////////////////////////////////////////

function intersect_edges(V::Points, edge1::Cell, edge2::Cell)
println("intersect_edges")
    err = 10e-8

    x1, y1, x2, y2 = vcat(map(c->V[c, :], edge1.nzind)...)
    x3, y3, x4, y4 = vcat(map(c->V[c, :], edge2.nzind)...)
    ret = Array{Tuple{Points, Float64}, 1}()

    v1 = [x2-x1, y2-y1];
    v2 = [x4-x3, y4-y3];
    v3 = [x3-x1, y3-y1];
    ang1 = dot(normalize(v1), normalize(v2))
    ang2 = dot(normalize(v1), normalize(v3))
    parallel = 1-err < abs(ang1) < 1+err
    colinear = parallel && (1-err < abs(ang2) < 1+err || -err < LinearAlgebra.norm(v3) < err)
    if colinear
        o = [x1 y1]
        v = [x2 y2] - o
        alpha = 1/dot(v,v')
        ps = [x3 y3; x4 y4]
        for i in 1:2
            a = alpha*dot(v',(reshape(ps[i, :], 1, 2)-o))
            if 0 < a < 1
                push!(ret, (ps[i:i, :], a))
            end
        end
    elseif !parallel
        denom = (v2[2])*(v1[1]) - (v2[1])*(v1[2])
        a = ((v2[1])*(-v3[2]) - (v2[2])*(-v3[1])) / denom
        b = ((v1[1])*(-v3[2]) - (v1[2])*(-v3[1])) / denom

        if -err < a < 1+err && -err <= b <= 1+err
            p = [(x1 + a*(x2-x1))  (y1 + a*(y2-y1))]
            push!(ret, (p, a))
        end
    end
    return ret
end

# //////////////////////////////////////////////////////////////////////////////

function spaceindex(model)::Array{Array{Int,1},1}
println("spaceindex")
	V,CV = model[1:2]
	dim = size(V,1)
	cellpoints = [ V[:,CV[k]]::Points for k=1:length(CV) ]
	#----------------------------------------------------------
	bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints]
	xboxdict = coordintervals(1,bboxes)
	yboxdict = coordintervals(2,bboxes)
	# xs,ys are IntervalTree type
	xs = IntervalTrees.IntervalMap{Float64, Array}()
	for (key,boxset) in xboxdict
		xs[tuple(key...)] = boxset
	end
	ys = IntervalTrees.IntervalMap{Float64, Array}()
	for (key,boxset) in yboxdict
		ys[tuple(key...)] = boxset
	end
	xcovers = boxcovering(bboxes, 1, xs)
	ycovers = boxcovering(bboxes, 2, ys)
	covers = [intersect(pair...) for pair in zip(xcovers,ycovers)]

	if dim == 3
		zboxdict = coordintervals(3,bboxes)
		zs = IntervalTrees.IntervalMap{Float64, Array}()
		for (key,boxset) in zboxdict
			zs[tuple(key...)] = boxset
		end
		zcovers = boxcovering(bboxes, 3, zs)
		covers = [intersect(pair...) for pair in zip(zcovers,covers)]
	end
	# remove each cell from its cover
	for k=1:length(covers)
		covers[k] = setdiff(covers[k],[k])
	end
	return covers
end

# //////////////////////////////////////////////////////////////////////////////

function cop2lar(cop::ChainOp)::Cells
println("cop2lar")
	[findnz(cop[k,:])[1] for k=1:size(cop,1)]
end

# //////////////////////////////////////////////////////////////////////////////
"""
    skel_merge(V1::Points, EV1::ChainOp, V2::Points, EV2::ChainOp)

Merge two **1-skeletons**
"""
function skel_merge(V1::Points, EV1::ChainOp, V2::Points, EV2::ChainOp)
println("skel_merge")
    V = [V1; V2]
    EV = blockdiag(EV1,EV2)
    return V, EV
end

"""
    skel_merge(V1::Points, EV1::ChainOp, FE1::ChainOp, V2::Points, EV2::ChainOp, FE2::ChainOp)

Merge two **2-skeletons**
"""
function skel_merge(V1::Points, EV1::ChainOp, FE1::ChainOp,
					V2::Points, EV2::ChainOp, FE2::ChainOp)
println("skel_merge")
    FE = blockdiag(FE1,FE2)
    V, EV = skel_merge(V1, EV1, V2, EV2)
    return V, EV, FE
end

# //////////////////////////////////////////////////////////////////////////////

function buildFV(EV::Cells, face::Cell)
println("buildFV")
    return buildFV(build_copEV(EV), face)
end

function buildFV(copEV::ChainOp, face::Cell)
println("buildFV")
    startv = -1
    nextv = 0
    edge = 0

    vs = Array{Int, 1}()

    while startv != nextv
        if startv < 0
            edge = face.nzind[1]
            startv = copEV[edge,:].nzind[face[edge] < 0 ? 2 : 1]
            push!(vs, startv)
        else
            edge = setdiff(intersect(face.nzind, copEV[:, nextv].nzind), edge)[1]
        end
        nextv = copEV[edge,:].nzind[face[edge] < 0 ? 1 : 2]
        push!(vs, nextv)

    end

    return vs[1:end-1]
end


# //////////////////////////////////////////////////////////////////////////////

function face_area(V::Points, EV::Cells, face::Cell)
println("face_area")
    return face_area(V, build_copEV(EV), face)
end

function face_area(V::Points, EV::ChainOp, face::Cell)
println("face_area")
    function triangle_area(triangle_points::Points)
println("triangle_area")
        ret = ones(3,3)
        ret[:, 1:2] = triangle_points
        return .5*det(ret)
    end

    area = 0

    fv = buildFV(EV, face)

    verts_num = length(fv)
    v1 = fv[1]

    for i in 2:(verts_num-1)

        v2 = fv[i]
        v3 = fv[i+1]

        area += triangle_area(V[[v1, v2, v3], :])
    end

    return area
end

# //////////////////////////////////////////////////////////////////////////////

function constrained_triangulation2D(V::Points, EV::Cells)
println("constrained_triangulation2D")
	triin = Triangulate.TriangulateIO()
	triin.pointlist = V
	triin.segmentlist = hcat(EV...)
	(triout, vorout) = Triangulate.triangulate("pQ", triin)
	trias = Array{Int64,1}[c[:] for c in eachcol(triout.trianglelist)]
	return trias
end

# //////////////////////////////////////////////////////////////////////////////

function triangulate2d(V::Points, EV::Cells)
println("triangulate2d")
   	 # data for Constrained Delaunay Triangulation (CDT)
   	 points = convert(Array{Float64,2}, V')
	 # points_map = Array{Int,1}(collect(1:1:size(points)[1]))
   	 # edges_list = convert(Array{Int,2}, hcat(EV...)')
   	 # edge_boundary = [true for k=1:size(edges_list,1)] ## dead code !!
	trias = constrained_triangulation2D(V::Points, EV::Cells)

 	#Triangle.constrained_triangulation(points,points_map,edges_list)
	innertriangles = Array{Int,1}[]
	for (u,v,w) in trias
		point = (points[u,:]+points[v,:]+points[w,:])./3
		copEV = lar2cop(EV)
		inner = point_in_face(point, points::Points, copEV::ChainOp)
		if inner
			push!(innertriangles,[u,v,w])
		end
	end
    return innertriangles
end

# //////////////////////////////////////////////////////////////////////////////

function planar_arrangement_2(V, copEV,bicon_comps, edge_map,
		sigma::Chain=spzeros(Int8, 0))
println("planar_arrangement_2")

    edges = sort(union(bicon_comps...))
    todel = sort(setdiff(collect(1:size(copEV,1)), edges))

    for i in reverse(todel)
        for row in edge_map

            filter!(x->x!=i, row)

            for j in 1:length(row)
                if row[j] > i
                    row[j] -= 1
                end
            end
        end
    end


	bicon_comps = biconnected_components(copEV)

	# component graph
	n, containment_graph, V, EVs, boundaries, shells, shell_bboxes = componentgraph(V,copEV,bicon_comps)

	copEV, FE = cell_merging(
	   	n, containment_graph, V, EVs, boundaries, shells, shell_bboxes)

	return V,copEV,FE
end

# //////////////////////////////////////////////////////////////////////////////

function biconnected_components(EV::ChainOp)
println("biconnected_components")

    ps = Array{Tuple{Int, Int, Int}, 1}()
    es = Array{Tuple{Int, Int}, 1}()
    todel = Array{Int, 1}()
    visited = Array{Int, 1}()
    bicon_comps = Array{Array{Int, 1}, 1}()
    hivtx = 1

    function an_edge(point) # TODO: fix bug
println("an_edge")
        # error? : BoundsError: attempt to access 0×0 SparseMatrix ...
        edges = setdiff(EV[:, point].nzind, todel)
        if length(edges) == 0
            edges = [false]
        end
        edges[1]
    end

    function get_head(edge, tail)
println("get_head")
        setdiff(EV[edge, :].nzind, [tail])[1]
    end

    function v_to_vi(v)
println("v_to_vi")
        i = findfirst(t->t[1]==v, ps)
        # seems findfirst changed from 0 to Nothing
        if typeof(i) == Nothing
            return false
        elseif i == 0
            return false
        else
            return ps[i][2]
        end
    end

    push!(ps, (1,1,1))
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
                for i in 1:size(EV,2)
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
                    edges = Array{Int, 1}()
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
    bicon_comps = sort(bicon_comps, lt=(x,y)->length(x)>length(y))
    return bicon_comps
end

# //////////////////////////////////////////////////////////////////////////////

function cleandecomposition(V, copEV, sigma, edge_map)
println("cleandecomposition")
    # Deletes edges outside sigma area
    todel = []
    new_edges = []
    map(i->new_edges=union(new_edges, edge_map[i]), sigma.nzind)
    ev = copEV[new_edges, :]
    for e in 1:copEV.m
        if !(e in new_edges)

            vidxs = copEV[e, :].nzind
            v1, v2 = map(i->V[vidxs[i], :], [1,2])
            centroid = .5*(v1 + v2)

            if ! point_in_face(centroid, V, ev)
                push!(todel, e)
            end
        end
    end

    for i in reverse(todel)
        for row in edge_map

            filter!(x->x!=i, row)

            for j in 1:length(row)
                if row[j] > i
                    row[j] -= 1
                end
            end
        end
    end

    V, copEV = delete_edges(todel, V, copEV)
	return V,copEV
end

# //////////////////////////////////////////////////////////////////////////////

function merge_vertices!(V::Points, EV::ChainOp, edge_map, err=1e-4)
println("merge_vertices")
    vertsnum = size(V, 1)
    edgenum = size(EV, 1)
    newverts = zeros(Int, vertsnum)
    # NearestNeighbors.KDTree constructor needs an explicit array of Float64
    V = Array{Float64,2}(V)
    kdtree = NearestNeighbors.KDTree(permutedims(V))

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
    edges = Array{Tuple{Int, Int}, 1}(undef, edgenum)
    oedges = Array{Tuple{Int, Int}, 1}(undef, edgenum)
    for ei in 1:edgenum
        v1, v2 = EV[ei, :].nzind
        edges[ei] = Tuple{Int, Int}(sort([newverts[v1], newverts[v2]]))
        oedges[ei] = Tuple{Int, Int}(sort([v1, v2]))
    end
    nedges = union(edges)
    nedges = filter(t->t[1]!=t[2], nedges)
    nedgenum = length(nedges)
    nEV = spzeros(Int8, nedgenum, size(nV, 1))
    # maps pairs of vertex indices to edge index
    etuple2idx = Dict{Tuple{Int, Int}, Int}()
    # builds `edge_map`
    for ei in 1:nedgenum
        nEV[ei, collect(nedges[ei])] .= 1
        etuple2idx[nedges[ei]] = ei
    end
    for i in 1:length(edge_map)
        row = edge_map[i]
        row = map(x->edges[x], row)
        row = filter(t->t[1]!=t[2], row)
        row = map(x->etuple2idx[x], row)
        edge_map[i] = row
    end
    # return new vertices and new edges
    return Points(nV), nEV
end

# //////////////////////////////////////////////////////////////////////////////

function frag_edge(V, EV::ChainOp, edge_idx::Int, bigPI)
println("frag_edge")
    alphas = Dict{Float64, Int}()
    edge = EV[edge_idx, :]
    verts = V[edge.nzind, :]
    for i in bigPI[edge_idx]
        if i != edge_idx
            intersection = intersect_edges(
            	V, edge, EV[i, :])
            for (point, alpha) in intersection
                verts = [verts; point]
                alphas[alpha] = size(verts, 1)
            end
        end
    end
    alphas[0.0], alphas[1.0] = [1, 2]
    alphas_keys = sort(collect(keys(alphas)))
    edge_num = length(alphas_keys)-1
    verts_num = size(verts, 1)
    ev = SparseArrays.spzeros(Int8, edge_num, verts_num)
    for i in 1:edge_num
        ev[i, alphas[alphas_keys[i]]] = 1
        ev[i, alphas[alphas_keys[i+1]]] = 1
    end
    return verts, ev
end

# //////////////////////////////////////////////////////////////////////////////

function planar_arrangement_1( V, copEV,
		sigma::Chain=spzeros(Int8, 0),
		return_edge_map::Bool=false,
		multiproc::Bool=false)
println("planar_arrangement_1")
	# data structures initialization
	edgenum = size(copEV, 1)
	edge_map = Array{Array{Int, 1}, 1}(undef,edgenum)
	rV = Points(zeros(0, 2))
	rEV = SparseArrays.spzeros(Int8, 0, 0)
	finalcells_num = 0

	# spaceindex computation
	model = (convert(Points,V'),cop2lar(copEV))
	bigPI = spaceindex(model)

    # multiprocessing of edge fragmentation
    if (multiproc == true)
        in_chan = Distributed.RemoteChannel(()->Channel{Int64}(0))
        out_chan = Distributed.RemoteChannel(()->Channel{Tuple}(0))
        ordered_dict = SortedDict{Int64,Tuple}()
        @async begin
            for i in 1:edgenum
                put!(in_chan,i)
            end
            for p in distributed.workers()
                put!(in_chan,-1)
            end
        end
        for p in distributed.workers()
            @async Base.remote_do(frag_edge_channel, p, in_chan, out_chan, V, copEV, bigPI)
        end
        for i in 1:edgenum
            frag_done_job = take!(out_chan)
            ordered_dict[frag_done_job[1]] = frag_done_job[2]
        end
        for (dkey, dval) in ordered_dict
            i = dkey
            v, ev = dval
            newedges_nums = map(x->x+finalcells_num, collect(1:size(ev, 1)))
            edge_map[i] = newedges_nums
            finalcells_num += size(ev, 1)
            rV, rEV = skel_merge(rV, rEV, v, ev)
        end
    else
        # sequential (iterative) processing of edge fragmentation
        for i in 1:edgenum
            v, ev = frag_edge(V, copEV, i, bigPI)
            newedges_nums = map(x->x+finalcells_num, collect(1:size(ev, 1)))
            edge_map[i] = newedges_nums
            finalcells_num += size(ev, 1)
            rV = convert(Points, rV)
            rV, rEV = skel_merge(rV, rEV, v, ev)
        end
    end
    # merging of close vertices and edges (2D congruence)
    V, copEV = rV, rEV
    V, copEV = merge_vertices!(V, copEV, edge_map)
	return V,copEV,sigma,edge_map
end

# //////////////////////////////////////////////////////////////////////////////

function permutationOrbits(perm::OrderedDict)
println("permutationOrbits")
	out = Array{Int64,1}[]
    while perm ≠ Dict()
        x = collect(keys(perm))[1]
        orbit = Int64[]
        while x in keys(perm)
            append!(orbit, perm[x])
			y,x = x,perm[x]
            delete!(perm,y)
		end
        append!(out, [ push!(orbit,orbit[1]) ]  )
	end
    return out
end

# //////////////////////////////////////////////////////////////////////////////

function faces2polygons(copEV,copFE)
println("faces2polygons")
	polygons = Array{Array{Int64,1},1}[]
	cycles = Array{Array{Array{Int64,1},1},1}[]
	for f=1:size(copFE,1)
		edges,signs = findnz(copFE[f,:])
		permutationMap = OrderedDict([ s>0 ? findnz(copEV[e,:])[1] : reverse(findnz(copEV[e,:])[1])
				for (e,s) in zip(edges,signs)])
		orbits = permutationOrbits(permutationMap)
		edgecycles = [[[ orbit[k], orbit[k+1] ] for k=1:length(orbit)-1]  for orbit in orbits]
		push!(polygons, [orbit[1:end-1] for orbit in orbits])
		push!(cycles, edgecycles)
	end
	return polygons,cycles
end

# //////////////////////////////////////////////////////////////////////////////

function triangulate2D(V::Points, cc::ChainComplex)::Array{Any, 1}
println("triangulate2D")
    copEV, copFE = cc
    triangulated_faces = Array{Any, 1}(undef, copFE.m)
    if size(V,2)==2
		V = [V zeros(size(V,1),1)]
	end

	polygons,edgecycles = faces2polygons(copEV, copFE) #new

    for f in 1:copFE.m
        edges_idxs = copFE[f, :].nzind
        edge_num = length(edges_idxs)
        edges = Array{Int,1}[] #zeros(Int, edge_num, 2)

		# fv = buildFV(copEV, copFE[f, :])
		fv = union(polygons[f]...)
        vs = V[fv, :]
		edges = union(edgecycles[f]...)
        edges = convert(Array{Int,2}, hcat(edges...)')

		# triangulated_faces[f] = Triangle.constrained_triangulation(
        # 	vs, fv, edges, fill(true, edge_num))
		v = convert(Points, vs'[1:2,:])
		vmap = Dict(zip(fv,1:length(fv))) # vertex map
		mapv = Dict(zip(1:length(fv),fv)) # inverse vertex map
		ev = [[vmap[e] for e in edges[k,:]] for k=1:size(edges,1)]
		trias = triangulate2d(v,ev)
		triangulated_faces[f] = [[mapv[v] for v in tria] for tria in trias]

        tV = V[:, 1:2]

        area = face_area(tV, copEV, copFE[f, :])
        if area < 0
            for i in 1:length(triangulated_faces[f])
                triangulated_faces[f][i] = triangulated_faces[f][i][end:-1:1]
            end
        end
    end

    return triangulated_faces
end

# //////////////////////////////////////////////////////////////////////////////

function FV2EVs(copEV::ChainOp, copFE::ChainOp)
println("FV2EVs")
	EV = [findnz(copEV[k,:])[1] for k=1:size(copEV,1)]
	FE = [findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
	EVs = [[EV[e] for e in fe] for fe in FE]
	return EVs
end

# //////////////////////////////////////////////////////////////////////////////

function planar_arrangement(
        V::Points,
        copEV::ChainOp,
        sigma::Chain=spzeros(Int8, 0),
        return_edge_map::Bool=false,
        multiproc::Bool=false)
println("planar_arrangement")

#planar_arrangement_1
	V,copEV,sigma,edge_map=planar_arrangement_1(V,copEV,sigma,return_edge_map,multiproc)

# cleandecomposition
	if sigma.n > 0
		V,copEV=cleandecomposition(V, copEV, sigma, edge_map)
	end

    bicon_comps = biconnected_components(copEV)
    # EV = cop2lar(copEV)
    # V,bicon_comps = biconnectedComponent((V,EV))

	if isempty(bicon_comps)
    	println("No biconnected components found.")
    	if (return_edge_map)
    	    return (nothing, nothing, nothing, nothing)
    	else
    	    return (nothing, nothing, nothing)
    	end
	end
#Planar_arrangement_2
	V,copEV,FE = planar_arrangement_2(V,copEV,bicon_comps,edge_map,sigma)
	if (return_edge_map)
	     return V, copEV, FE, edge_map
	else
	     return V, copEV, FE
	end
end

# //////////////////////////////////////////////////////////////////////////////

function characteristicMatrix( FV::Cells )::ChainOp
println("characteristicMatrix")
	I,J,V = Int64[],Int64[],Int8[]
	for f=1:length(FV)
		for k in FV[f]
		push!(I,f)
		push!(J,k)
		push!(V,1)
		end
	end
	M_2 = sparse(I,J,V)
	return M_2
end

# //////////////////////////////////////////////////////////////////////////////

function boundary_1( EV::Cells )::ChainOp
println("boundary_1")
	out = characteristicMatrix(EV)'
	for e = 1:length(EV)
		out[EV[e][1],e] = -1
	end
	return out
end

coboundary_0(EV::Cells) = convert(ChainOp,LinearAlgebra.transpose(boundary_1(EV::Cells)))

# //////////////////////////////////////////////////////////////////////////////

function arrange2D(V,EV)
println("arrange2D")
	cop_EV = coboundary_0(EV::Cells)
	cop_EW = convert(ChainOp, cop_EV)
	W = convert(Points,V')
	V, copEV, copFE = planar_arrangement(W::Points, cop_EW::ChainOp)
	EVs = FV2EVs(copEV, copFE) # polygonal face fragments

	triangulated_faces = triangulate2D(V, [copEV, copFE])
	FVs = convert(Array{Cells}, triangulated_faces)
	V = convert(Points,V')
	return V,FVs,EVs, copEV, copFE
end

# //////////////////////////////////////////////////////////////////////////////

function removedups(obj)::Cells
println("removedups")
   # initializations
   hulls = ToLAR(obj).childs[1].facets
   dict = ToLAR(obj).childs[1].db
   inverted_dict = Dict{valtype(dict), Vector{keytype(dict)}}()
   [push!(get!(() -> valtype(inverted_dict)[], inverted_dict, v), k) for (k, v) in dict]  
   DB = []  # convert arrays of indices to arrays of points
   for hull in hulls
      points = []
      [ append!(points, inverted_dict[k]) for k in hull ]
      push!(DB, Set(points))
   end 
   DB = Set(DB) # eliminate duplicates
   faces = [[dict[point] for point in set] for set in DB]
   faces = sort!(AA(sort!)(faces))
end

function LAR(obj)
println("LAR")
   V = ToLAR(obj).childs[1].points
   CV = ToLAR(obj).childs[1].hulls
   facets = ToLAR(obj).childs[1].facets
   FV = removedups(obj)
   FF = CSC(FV) * CSC(FV)'
   edges = filter(x->x[1]<x[2] && FF[x...]==2,collect(zip(findnz(FF)[1:2]...)))
   EV = sort!(collect(Set([FV[i] ∩ FV[j] for (i,j) in edges])))
   out = Plasm.Lar(hcat(V...), Dict(:CV=>CV, :FV=>FV, :EV=>EV))
end



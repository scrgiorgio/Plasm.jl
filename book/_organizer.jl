#///// 7 Hierarchical Structures   /////////////////////


#///// 7.1 Assembly of Components   /////////////////////
#///// 7.1.1 Recursive Bill of Materials (BOM)  /////////////////////
#///// 7.1.2 Model (object) Transform   /////////////////////
#///// 7.1.3 Viewing (camera) Transform   /////////////////////


#///// 7.2 Traversal Algorithms  /////////////////////
#///// 7.2.1 Structure/Component Definition   /////////////////////
#///// 7.2.2 Mapping local to global coordinates  /////////////////////
#///// 7.2.3 Abstract Syntax Tree (AST) Linearize   /////////////////////


# //////////////////////////////////////////////////////////////////////////////

	using DataStructures
	using IntervalTrees
   using Triangulate
   using Triangle
	using SparseArrays
	using LinearAlgebra
	using NearestNeighbors
   using ViewerGL;   GL = ViewerGL

	const Points = Matrix
	const Cells = Vector{Vector{Int}}
	const Cell = SparseVector{Int8, Int}
	const Chain = SparseVector{Int8,Int}
	const ChainOp = SparseMatrixCSC{Int8,Int}
	const ChainComplex = Vector{ChainOp}
	const LARmodel = Tuple{Points, Vector{Cells}}
	const LAR = Union{ Tuple{Points, Cells}, Tuple{Points, Cells, Cells}, Tuple{Points, Cells, Cells, Cells} }
	Verbose = false




#///// 8 Intersection of Input Cells   /////////////////////


#///// 8.1.1 Outline of face splitting   /////////////////////

#///// 8.1.2 Computing a 2D face arrangement  /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Test for containment of checkpoint into triangle `f`; used in TGW algorithm """
function interior_to_f(triangle,f,V,FV,EV,FE)::Bool
	v1,v2,v3 = triangle
	u = V[:,v2]-V[:,v1]
	v = V[:,v3]-V[:,v1]
	w = cross(u,v)
	T = [1. 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; T[1:3,4] = -V[:,v1]
	R = [1. 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; R[1:3,1:3] = [u v w]; R = R'
	mapping = R * T

	trianglepoints = [[V[:,v1] V[:,v2] V[:,v3]]; ones(3)']
	pts = mapping * trianglepoints
	points2D = [pts[r,:] for r = 1:size(pts,1)
				if !(all(pts[r,:].==0) || all(pts[r,:].==1) )]
	p2D = hcat(points2D...)'
	checkpoint = 0.4995 .* p2D[:,1] + 0.4995 .* p2D[:,2] + 0.001 .* p2D[:,3]

	cellpoints = [V[:,FV[f]]; ones(length(FV[f]))' ]
	points = mapping * cellpoints
	verts2D = [points[r,:] for r = 1:size(points,1)
				if !(all(points[r,:].==0) || all(points[r,:].==1) )]
	#P2D = hcat(verts2D...)'
	P2D = hcat(verts2D...)
   
	vdict = Dict(collect(zip(FV[f], 1:length(FV[f]))))
	celledges = [[vdict[v] for v in EV[e]] for e in FE[f]]

   inner = point_in_face(checkpoint, P2D::Points, lar2cop(celledges)::ChainOp)
end


#///// 8.2 Accelerating the search for intersections   /////////////////////

#///// 8.2.1 Boxing input (d−1)-cells in dD   /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Test for relative containment of two boxes; used with 2D graphs """
function bbox_contains(container, contained)::Bool
    b1_min, b1_max = container
    b2_min, b2_max = contained
    all(map((i,j,k,l) -> i<=j<=k<=l, b1_min, b2_min, b2_max, b1_max))
end

""" Multidimensional constructor of bounding box of a vertex subset """
function bbox(vertices::Points)
    minimum = mapslices(x->min(x...), vertices, dims=1)
    maximum = mapslices(x->max(x...), vertices, dims=1)
    minimum, maximum
end

#///// 8.2.2 Using d interval trees   /////////////////////
# //////////////////////////////////////////////////////////////////////////////
using IntervalTrees

""" Compute the boundingbox for multidimensional points """
function boundingbox(vertices::Points)
   minimum = mapslices(x->min(x...), vertices, dims=2)
   maximum = mapslices(x->max(x...), vertices, dims=2)
   return minimum, maximum
end

""" Make dictionary of boxes for IntervalTrees construction """
function coordintervals(coord, bboxes)	
	boxdict = OrderedDict{Vector{Float64}, Vector{Int64}}()
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

""" Compute the `tree` 1D `covers` for each bounding box `extent` """
function boxcovering(bboxes, index, tree)
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

""" Main to build the multidimensional `spaceindex(model)` for any d-cell type """
function spaceindex(model)::Vector{Vector{Int}}
	V,CV = model[1:2]
	dim = size(V,1)
	cellpoints = [ V[:,CV[k]]::Points for k=1:length(CV) ]
	bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints]
	boxdicts = [coordintervals(k,bboxes) for k=1:dim]
	trees = [IntervalTrees.IntervalMap{Float64, Array}() for k=1:dim]
	covers = []
	for k = 1:dim
		for (key,boxset) in boxdicts[k]
			trees[k][tuple(key...)] = boxset
		end
		push!(covers, boxcovering(bboxes, k, trees[k]))
	end
	covers = [reduce(intersect,[covers[h][k] for h=1:dim]) for k=1:length(CV)]
end


#///// 8.2.3 Intersect query results   /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Arrangement of 2D space starting from 1-skeleton in 2D """
function planar_arrangement(
        V::Points,
        copEV::ChainOp,
        sigma::Chain=spzeros(Int8, 0),
        return_edge_map::Bool=false)
# `sigma::Chain`: if specified, `planar_arrangement` will delete from the output every edge and face outside this cell. Defaults to an empty cell.
# `return_edge_map::Bool`: makes the function return also an `edge_map` which maps the edges of the imput to the one of the output. Defaults to `false`.

# task planar_arrangement_1
	V,copEV,sigma,edge_map = planar_arrangement_1(V,copEV,sigma,return_edge_map)

# task cleandecomposition
	if sigma.n > 0
		V,copEV = cleandecomposition(V, copEV, sigma, edge_map)
	end

   #bicon_comps = biconnectedComponent(cop2lar(copEV))
   bicon_comps = biconnected_components(copEV)
#@show bicon_comps;
	if isempty(bicon_comps)
    	println("No biconnected components found.")
    	if (return_edge_map)
    	    return (nothing, nothing, nothing, nothing)
    	else
    	    return (nothing, nothing, nothing)
    	end
	end
# task Planar_arrangement_2
	V,copEV,FE = planar_arrangement_2(V,copEV,bicon_comps,edge_map,sigma)
	if (return_edge_map)
	     return V, copEV, FE, edge_map
	else
	     return V, copEV, FE
	end
end

# //////////////////////////////////////////////////////////////////////////////

""" Build the 2D 1-skeleton of arranged face `sigma` """
function planar_arrangement_1( V, copEV,
		sigma::Chain=spzeros(Int8, 0),
		return_edge_map::Bool=false )
	# data structures initialization
	edgenum = size(copEV, 1)
	edge_map = Vector{Vector{Int}}(undef,edgenum)
	rV = Points(zeros(0, 2))
	rEV = SparseArrays.spzeros(Int8, 0, 0)
	finalcells_num = 0

	# spaceindex computation
	model = (convert(Points,V'), cop2lar(copEV))
	bigPI = spaceindex(model::LAR)

    # sequential (iterative) processing of edge fragmentation
   for i in 1:edgenum
      v, ev = frag_edge(V, copEV, i, bigPI)
      newedges_nums = map(x->x+finalcells_num, collect(1:size(ev, 1)))
      edge_map[i] = newedges_nums
      finalcells_num += size(ev, 1)
      rV = convert(Points, rV)
      rV, rEV = skel_merge(rV, rEV, v, ev)
   end
    # merging of close vertices and edges (2D congruence)
    V, copEV = rV, rEV
    V, copEV = merge_vertices!(V, copEV, edge_map)
	return V,copEV,sigma,edge_map
end


# //////////////////////////////////////////////////////////////////////////////

""" Final task of 2D `planar_arrangement` """
function planar_arrangement_2(V, copEV,bicon_comps, edge_map,
		   sigma::Chain=spzeros(Int8, 0))

   #edges = sort(union(bicon_comps...))
   edges = union(bicon_comps...)
   #todel = sort(setdiff(collect(1:size(copEV,1)), edges))
   todel = setdiff(collect(1:size(copEV,1)), edges)

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
	#bicon_comps = biconnectedComponent(cop2lar(copEV))
	
	# component graph
	n, containment_graph, V, EVs, boundaries, shells, shell_bboxes = 
	      componentgraph(V, copEV, bicon_comps)
	copEV, FE = cell_merging(n, containment_graph, 
	      V, EVs, boundaries, shells, shell_bboxes)

	return V,copEV,FE
end

#///// 8.3 2-cell Splitting   /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Merge two 1-chains; used by `frag_face` """
function skel_merge(V1::Points, EV1::ChainOp, V2::Points, EV2::ChainOp)
   V = [V1; V2]
   EV = blockdiag(EV1,EV2)
   return V, EV
end

# //////////////////////////////////////////////////////////////////////////////

""" taglia con z=0 la i-ma faccia di sp_idx[sigma] """
# eliminando il `for` esegue tutto l'arrangement 
# ma non esegue lo splitting di (ogni) `face` 
# quindi riporta i dati di input alla pipeline !!

function face_int(V::Points, EV::ChainOp, face::Cell)
    vs = buildFV(EV, face)     # EV::ChainOp, face::Cell
#@show vs;
    retV = Points(undef, 0, 3)
    err = 10e-8
    visited_verts = []
   for i in 1:length(vs)
       o = V[vs[i],:]
       j = i < length(vs) ? i+1 : 1
       d = V[vs[j],:] - o    # vertex in local coordinates
       if !(-err < d[3] < err)
           alpha = -o[3] / d[3]
           if -err <= alpha <= 1+err
               p = o + alpha*d
               if -err < alpha < err || 1-err < alpha < 1+err
                   if !(vin(p, visited_verts))
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
#@show retV, retEV;
    retV, retEV
# ESEMPIO:
# vs = [19, 20, 24, 23]
# (retV, retEV) = (Matrix{Any}(undef, 0, 3), sparse(Int64[], Int64[], Int8[], 0, 0))
# vs = [17, 19, 23, 21]
# (retV, retEV) = (Any[0.01580735103001162 0.28058426264671665 0.0; -0.054613822013010384 0.32096378672306647 0.0], sparse([1, 1], [
end

# //////////////////////////////////////////////////////////////////////////////

""" Splitting of `sigma` face against faces in `sp_idx[sigma]` """
function frag_face(V, EV::ChainOp, FE::ChainOp, sp_idx, sigma)

   vs_num = size(V, 1)

# 2D transformation of `sigma` face
   sigmavs = (abs.(FE[sigma:sigma,:]) * abs.(EV))[1,:].nzind # sigma vertex indices
   sV = V[sigmavs, :]
   sEV = EV[FE[sigma, :].nzind, sigmavs]
   M = submanifold_mapping(sV)
   tV = ([V ones(vs_num)]*M)[:, 1:3]  # folle convertire *tutti* i vertici
   sV = tV[sigmavs, :]
   # `sigma` face intersection with faces in `sp_idx[sigma]`, i.e., in `bigpi`
   for i in sp_idx[sigma]
       tmpV, tmpEV = face_int(tV, EV, FE[i, :])
       sV, sEV = skel_merge(sV, sEV, tmpV, tmpEV)
   end
   # computation of 2D arrangement of sigma face
   sV = sV[:, 1:2]
   nV, nEV, nFE = planar_arrangement(sV, sEV, sparsevec(ones(Int8, length(sigmavs))))
   nvsize = size(nV, 1)
   # return each 2D complex in 3D
   nV = [nV zeros(nvsize) ones(nvsize)] * inv(M)[:, 1:3] 
   return nV, nEV, nFE
end

## ////////////////////////////////////////////////////////////////////////////// incomleta
#
#"""
#	sigma_intersect(V, EV, FE, sigma, Q, bigpi)
#
#Intersect `sigma` plane `z=0` with edges in transformed `bigpi`.
#"""
#function sigma_intersect(V, EV, FE, sigma, Q, bigpi)
#	sigma_edges = EV[FE[sigma]]
#	sigma_verts = union(sigma_edges...)
#	sigma_vdict =  Dict(zip(sigma_verts, 1:length(sigma_verts)))
#	sigma_lines = [[sigma_vdict[v] for v in edge] for edge in sigma_edges]
#	# sigma face on z=0 plane
#	S = V[:,sigma_verts] # sigma vertices
#	Z = (Q * [S; ones(1,size(S,2))])[1:3,:] # sigma mapped in z=0
#	#Plasm.view(Z, sigma_lines)
#	# initialization of line storage (to be intersected)
#	linestore = [[Z[:,v1] Z[:,v2]] for (v1,v2) in sigma_lines]
#	linenum = length(sigma_lines)
#	# reindexing of ∏(σ) chain
#	bigpi_edges = [EV[fe] for fe in FE[bigpi]]
#	bigpi_verts = union(union(bigpi_edges...)...)
#	bigpi_vdict =  Dict(zip(bigpi_verts, 1:length(bigpi_verts)))
#	bigpi_lines = [[sort([bigpi_vdict[v] for v in edge]) for edge in faceedges]
#		for faceedges in bigpi_edges]
#	# bigpi trasformed in 3D according to Q mapping  ==> BUG ??!! :  NO
#	P = V[:,bigpi_verts] # bigpi vertices
#	W = (Q * [P; ones(1,size(P,2))])[1:3,:] # bigpi mapped by Q
#	#Plasm.view(W, union(bigpi_lines...))
#	# filter on bigpi_lines that do not cross z=0
#	filtered_edges = [[[v1,v2] for (v1,v2) in face_lines
#		if (sign(W[3,v1]) * sign(W[3,v2])) <= 0]
#			for face_lines in bigpi_lines]
#	filtered_edges = [edge for edge in filtered_edges if edge!=[]]
#	#Plasm.view(W, union(filtered_edges...))
#	#Plasm.view(Plasm.numbering()((W,[[[k] for k=1:size(W,2)],union(filtered_edges...)])))
#	# computation of ordered z=0 points by face
#	facepoints = []
#	for face in filtered_edges
#		points = Array{Float64,1}[]
#		for (v1,v2) in face
#			#v1[3] + (v2[3]-v1[3])*t = 0 # compute parameter of intersection point with z=0
#			if (W[3,v2] - W[3,v1]) != 0
#				t = -W[3,v1] / (W[3,v2] - W[3,v1])  # z1 + (z2-z1)*t = 0
#				point = W[:,v1] + (W[:,v2] - W[:,v1]) * t
#				push!(points,point)
#			end
#		end
#		if length(Set(points)) > 1
#			push!(facepoints,Set(points))
#		end
#		# TODO: linearly order facepoints[face]
#	end
#	#Plasm.view(Lar.Struct([ (hcat(collect(edge)...), [[1,2]]) for edge in facepoints ]))
#	# compute z_lines
#	c = collect
#	z0_lines = [[[c(ps)[k] c(ps)[k+1]] for k=1:2:(length(ps)-1)] for ps in facepoints]
#
#	# intersecting lines upload in linestore
#	linestore = [[Z[:,v1] Z[:,v2]] for (v1,v2) in sigma_lines]
#	if z0_lines != []
#		linestore = append!(linestore, union(z0_lines...))
#	end
#	return linestore,linenum
#end
#
## ////////////////////////////////////////////////////////////////////////////// incomleta
#
#"""
#	computeparams(linestore,linenum)
#
#Intersect all (non-simmetric) pairs of `line` in `linestore`.
#"""
#function computeparams(linestore,linenum)
#	m = length(linestore)
#	params = [[] for i=1:m]
#	for h=1:(m-1)
#		(x1,y1,x2,y2) = reshape(linestore[h][1:2,:],1,4)
#		for k=(h+1):m
#			(x3,y3,x4,y4) = reshape(linestore[k][1:2,:],1,4)
#			# intersect lines h,k
#			# http://www.cs.swan.ac.uk/~cssimon/line_intersection.html
#			det = (x4-x3)*(y1-y2)-(x1-x2)*(y4-y3)
#			α,β = 99,99
#			if det != 0.0
#				a = 1/det
#				b = [y1-y2 x2-x1; y3-y4 x4-x3]  # x1-x2 => x2-x1 bug in the source link !!
#				c = [x1-x3; y1-y3]
#				(β,α) = a * b * c
#			else
#				if (y1==y2) == (y3==y4) # segments collinear
#					α = -x1/(x2-x1)
#					β = -x3/(x4-x3)
#				elseif (x1==x2) == (x3==x4)
#					α = -y1/(y2-y1)
#					β = -y3/(y4-y3)
#				else
#					 # segments parallel: no intersection
#					 (β,α) = 999
#				end
#			end
#			if 0<α<1 && 0<β<1
#				push!(params[h], α)
#				push!(params[k], β)
#			end
#		end
#	end
#	out = Array{Float64,1}[]
#	for line in params
#		push!(line, 0.0)
#		push!(line, 1.0)
#		line = sort(collect(Set(line)))
#		push!(out, line)
#	end
#	return out
#end
#
## ////////////////////////////////////////////////////////////////////////////// incomleta
#
#""" Splitting of `sigma` face against faces in `bigpi`, i.e., `sp_idx[sigma]` """
#function fragface(V, cop_EV, cop_FE, sp_idx, sigma)
#	# format conversion of function parameters
#	EV = [findnz(cop_EV[k,:])[1] for k=1:size(cop_EV,1)]
#	FE = [findnz(cop_FE[k,:])[1] for k=1:size(cop_FE,1)]
#   FV = [union(cat([EV[e] for e in FE[f]],dims=2)...) for f=1:length(FE)]
#	V = convert(Points, V')
#
#	#Q = face_mapping(V, FV, sigma)
#	Q = facemapping(V, FV, sigma)
#	bigpi = sp_idx[sigma]
#	linestore, linenum = sigma_intersect(V, EV, FE, sigma, Q, bigpi)
#@show linestore, linenum;	
#	lineparams = computeparams(linestore,linenum)
#	pairs = collect(zip(lineparams,linestore))
#@show pairs;
#
#	linepoints = Vector{Vector{Float64}}[]
#	number = 0
#	for (k,(params,line)) in enumerate(pairs)
#		v1,v2 = line[1:2,1],line[1:2,2]
#		if params==[0.0, 1.0]
#			points = [v1,v2]
#			if k<=linenum number += 1 end
#		elseif length(params)>2
#			points = [v1+α*(v2-v1) for α ∈ params]
#			if k<=linenum number += length(params)-1 end
#		end
#		push!(linepoints,points)
#	end
#	verts = []
#	edges = Vector{Float64}[]
#	offset = 0
#	for (h,pointline) in enumerate(linepoints)
#		offset
#		append!(verts, pointline)
#		append!(edges, [[k+offset,k+1+offset] for k=1:length(pointline)-1])
#		offset = length(verts)
#	end
#	verts = hcat(verts...)
#	verts = convert(Points, verts')
#
#	copEV = Lar.coboundary_0(convert(Cells,edges))
#	copFE = coboundary_1(verts,faces,edges)
#	return verts, copEV, copFE
#end

#///// 8.3.1 Mapping2-cells to z=0  /////////////////////
# //////////////////////////////////////////////////////////////////////////////KO incomleta
""" Compute the map from vs (at least three) to z=0 """
# REMARK: will not work when vs[:,1:3] are aligned !!!!  TODO: fix 
function submanifold_mapping(vs)
   #@show vs;
    u1 = vs[2,:] - vs[1,:]
    u2 = vs[3,:] - vs[1,:]
    u3 = LinearAlgebra.cross(u1, u2)
    T = Matrix{Float64}(LinearAlgebra.I, 4, 4)
    T[4, 1:3] = - vs[1,:]
    M = Matrix{Float64}(LinearAlgebra.I, 4, 4)
    M[1:3, 1:3] = [u1 u2 u3]
    return T * M
end

# //////////////////////////////////////////////////////////////////////////////KO incomleta

#""" aaaaaaaaa """
#function facemapping(V, FV, sigma; approx=10e-8 )
#	vs = FV[sigma]; i = 1
#	# compute affinely independent triple
#	n = length(vs)
#	succ(i) = i % n + 1
#	a = normalize(V[:,vs[succ(i)]] - V[:,vs[i]])
#	b = normalize(V[:,vs[succ(succ(i))]] - V[:,vs[succ(i)]])
#	c = cross(a,b)
#	while (-approx < det([a b c]) < approx)
#		i += 1
#		a = normalize(V[:,vs[succ(i)]] - V[:,vs[i]])
#		b = normalize(V[:,vs[succ(succ(i))]] - V[:,vs[succ(i)]])
#		c = cross(a,b)
#	end
#	# sigma translation
#	T = Lar.t(-V[:,vs[succ(i)]]...)
#	W = T * [V[:,vs] ; ones(1,length(vs))]
#	t_vs = W[1:3,:]
#	# translated sigma rotation
#	r1,r2,r3 = t_vs[:,i], t_vs[:,succ(succ(i))], cross(t_vs[:,i], t_vs[:,succ(succ(i))])
#	eye(n) = LinearAlgebra.Matrix{Float64}(I,n,n)
#	R = eye(4)
#	R[1:3,1:3] = inv([r1 r2 r3])
#	mapping = R * T
#	return mapping
#end


#///// 8.3.2 Intersecting candidate planes with z=0 /////////////////////

#///// 8.3.3 Mutually splitting 1-cells in 2D   /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Intersection of two edges in 2D space (Maria Teresa?) """
function intersect_edges(V::Points, edge1::Cell, edge2::Cell)
    err = 10e-8
#@show V, findnz(edge1), findnz(edge2);
    x1, y1, x2, y2 = vcat(map(c->V[c, :], edge1.nzind)...)
    x3, y3, x4, y4 = vcat(map(c->V[c, :], edge2.nzind)...)
    ret = Array{Tuple{Points, Float64}, 1}()

    v1 = [x2-x1, y2-y1];
    v2 = [x4-x3, y4-y3];
    v3 = [x3-x1, y3-y1];
    ang1 = dot(normalize(v1), normalize(v2))
    ang2 = dot(normalize(v1), normalize(v3))
    parallel = 1-err < abs(ang1) < 1+err
    colinear = parallel && (1-err < abs(ang2) < 1+err || -err < norm(v3) < err)
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
	
""" Fragment in 2D the `edge_idx` edge """
function frag_edge(V, EV::ChainOp, edge_idx::Int, bigPI)
   #@show V, findnz(EV), edge_idx, bigPI[edge_idx];
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

#///// 8.4 Local topology computation   /////////////////////

# //////////////////////////////////////////////////////////////////////////////

""" Find external cycles (area < 0) in chain complex (FE already computed) """
function get_external_cycle(V::Points, EV::ChainOp, FE::ChainOp)
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

""" 2D partial input (`function edge_angle`) for multidimensional TGW """
function minimal_2cycles(V::Points, EV::ChainOp)

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
    VE = convert(ChainOp, EV')
    EF = minimal_cycles(edge_angle)(V, VE) #  call the 2D TGW implementation

    return convert(ChainOp, EF')
end

#///// 8.4.1 Extract maximal biconnected graphs  /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Prune the 2D (graph) of local complex in `planar_arrangement_2`  """
function componentgraph(V, copEV, bicon_comps)
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

#///// 8.4.2 Topological gift wrapping in 2D   /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Main for the Furiani's dimension-independent (2D-3D) TGW; used only in 2D """
function minimal_cycles(angles_fn::Function)

    """ Actual computational engine for 2D TGW algorithm """
    function _minimal_cycles(V::Points, ld_bounds::ChainOp)  #  EV)
    
        # sembra scritto appositamente per essere illegibile !!
        # posso riscriverlo in modo più umano.
        # molto dipende dall'uso di matrici sparse, inutili, in quanto i
        # main loops non sono algebrici.  Giorgio: Fammi sapere.
        # Si può adattare a 2D l'implementazione 3D (mia) molto più semplice

        lld_cellsnum, ld_cellsnum = size(ld_bounds)
        count_marks = zeros(Int64, ld_cellsnum)
        dir_marks = zeros(Int64, ld_cellsnum)
        d_bounds = spzeros(Int64, ld_cellsnum, 0)

        angles = Array{Array{Int64, 1}, 1}(undef, lld_cellsnum)

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
            sort!(as, lt=(a,b)->a[2]<b[2])
            as = map(a->a[1], as)
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
            for k=1:length(count_marks)
            	if count_marks[k]>2  error("TGW is looping") end
            end

            as[ne]
        end

        while (sigma = get_seed_cell()) > 0
            c_ld = spzeros(Int8, ld_cellsnum)
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

#///// 8.4.3 Set of 2D local chain complexes   /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Delete graph edges outside `sigma`s areas; used by `cleandecomposition` """
function delete_edges(todel, V::Points, EV::ChainOp)
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

""" Remove external incomplete (open) parts from 2D local arrangements  """
function cleandecomposition(V, copEV, sigma, edge_map)
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
	
""" Probably not useful; used in `componentgraph`: Quadratic in 2D shells !! 
Implementation of DAG of relative containment of complexes (ACM TSAS) """
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
	
""" Implementation of DAG of relative containment of complexes (ACM TSAS) """
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
	
""" Implementation of DAG of relative containment of complexes (ACM TSAS) """
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

#///// 0.0.0 Utilities   /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Collect the edges into different subsets for each face in `copFE` """
function FV2EVs(copEV::ChainOp, copFE::ChainOp)
	EV = [findnz(copEV[k,:])[1] for k=1:size(copEV,1)]
	FE = [findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
	EVs = [[EV[e] for e in fe] for fe in FE]
	return EVs
end

# //////////////////////////////////////////////////////////////////////////////

""" Return the Lar `Cells` (Vector{Vector{Int64}}) datatype """
function cop2lar(cop::ChainOp)::Cells
	[findnz(cop[k,:])[1] for k=1:size(cop,1)]::Cells
end

# //////////////////////////////////////////////////////////////////////////////

""" Component of TGW in 3D (Pao);  return an ordered `fan` of 2-cells """
function ord(hinge::Int, bd1::AbstractSparseVector{Int,Int}, V::Points,
             FV::Cells, EV::Cells, FE::Cells)
   
	cells = SparseArrays.findnz(bd1)[1]
	triangles = []

	function area(v1,v2,v3)
		u = V[:,v2]-V[:,v1]
		v = V[:,v3]-V[:,v1]
		out = norm(cross(u,v)) # actually, to be divided by two
		return out
	end

	for f in cells
		v1,v2 = EV[hinge]
		index = findfirst(v -> (area(v1,v2,v)≠0), FV[f])
		v3 = FV[f][index]

		# test if [v1,v2,v3] interior to f
		while true
			if interior_to_f([v1, v2, v3],f,V,FV,EV,FE)
				push!(triangles, [v1,v2,v3])
				break
			else
				index = findnext(v -> (area(v1,v2,v)≠0), FV[f], index+1)
				v3 = FV[f][index]
			end
		end
	end
	order = ordering(triangles,V)
	return [cells[index] for index in order]
end

# //////////////////////////////////////////////////////////////////////////////

""" Alternate method for buildFV function """
function buildFV(EV::Cells, face::Cell)
    #return buildFV(build_copEV(EV), face)
    return buildFV(lar2cop(EV), face)
end

""" Alternate method for buildFV function """
# REMARK:  what method is used where ??
function buildFV(copEV::ChainOp, face::Cell)
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

""" Alternate method for `face_area` function """
function face_area(V::Points, EV::Cells, face::Cell)
    return face_area(V, build_copEV(EV), face)
end

""" Basic algorithm for polygon area by accumulating signed triangle areas """
function face_area(V::Points, EV::ChainOp, face::Cell)
    
    """ Compute signed triangle area from vertex points """
    function triangle_area(triangle_points::Points)
        ret = ones(3,3)
        ret[:, 1:2] = triangle_points
        return .5*det(ret)
    end

    area = 0

    fv = buildFV(EV, face) # EV::Cells

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

""" Predicate to check membership of `vertex` in `vertices_set` array"""
function vin(vertex, vertices_set)::Bool
    for v in vertices_set
        if vequals(vertex, v)
            return true
        end
    end
    return false
end

# //////////////////////////////////////////////////////////////////////////////

""" Predicate to check equality of two vertices (only used above) """
function vequals(v1, v2)
    err = 10e-8
    return length(v1) == length(v2) && all(map((x1, x2)->-err < x1-x2 < err, v1, v2))
end

# //////////////////////////////////////////////////////////////////////////////

""" conversion from Lar `Cells` to Lar `ChainOp` """
function lar2cop( CV::Cells )::ChainOp
	I = vcat( [ [k for h in CV[k]] for k=1:length(CV) ]...)	
	J = vcat( CV...)	
	X = Int8[1 for k=1:length(I)] 	
	CV = SparseArrays.sparse(I,J,X)	
end

# //////////////////////////////////////////////////////////////////////////////

""" conversion from Lar `ChainOp` to Lar `Cells` """
function cop2lar(CV::ChainOp)::Cells
	CV = [findnz(CV[k,:])[1] for k=1:size(CV,1)]
end

# //////////////////////////////////////////////////////////////////////////////

""" conversion from 3D `ChainComplex` to `Vector` of 1-chains in 3D """
function FV2EVs(copEV::ChainOp, copFE::ChainOp)
	EV = [findnz(copEV[k,:])[1] for k=1:size(copEV,1)]
	FE = [findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
	EVs = [[EV[e] for e in fe] for fe in FE]
	return EVs
end

#///// 9 Congruence of Local Chains   /////////////////////
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
   W = convert(Points, LinearAlgebra.transpose(V))
   kdtree = KDTree(W)
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
   V[30,:] = V[26,:] #<<<<<<<<<<<<<<<<<< TEST

   # translate edges to take congruence into account
   edges = Array{Tuple{Int, Int}, 1}(undef, edgenum)
   oedges = Array{Tuple{Int, Int}, 1}(undef, edgenum)
   for ei in 1:edgenum
       v1, v2 = EV[ei, :].nzind
       edges[ei] = Tuple{Int, Int}(sort([newverts[v1], newverts[v2]]))
       oedges[ei] = Tuple{Int, Int}(sort([v1, v2]))
   end
   nedges = union(edges)
   # remove edges of zero length
   nedges = filter(t->t[1]!=t[2], nedges)
   nedgenum = length(nedges)
   nEV = spzeros(Int8, nedgenum, size(nV, 1))

   etuple2idx = Dict{Tuple{Int, Int}, Int}()
   for ei in 1:nedgenum
   	begin
       	nEV[ei, collect(nedges[ei])] .= 1
       	nEV
       end
       etuple2idx[nedges[ei]] = ei
   end
   for e in 1:nedgenum
   	v1,v2 = findnz(nEV[e,:])[1]
   	nEV[e,v1] = -1; nEV[e,v2] = 1
   end

   # compute new faces to take congruence into account
   faces = [[
       map(x->newverts[x], FE[fi, ei] > 0 ? oedges[ei] : reverse(oedges[ei]))
       for ei in FE[fi, :].nzind
   ] for fi in 1:facenum]


   visited = []
   function filter_fn(face)

       verts = []
       map(e->verts = union(verts, collect(e)), face)
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
           ei = etuple2idx[Tuple{Int, Int}(sort(collect(edge)))]
           nFE[fi, ei] = sign(edge[2] - edge[1])
       end
   end

   return Points(nV), nEV, nFE
end


#///// 9.1 Quotient topology  167  /////////////////////
#///// 9.1.1 Nearness relation167  /////////////////////
#///// 9.1.2 Topological congruence  /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Search maximal biconnected components in a graph (no pao) """
# Hopcroft-Tarjan algorithm. Depth-First-Visit recursive algorithm.
function DFV_visit( VV::Cells, out::Array, count::Int, visited::Array, parent::Array, d::Array, low::Array, stack::Array, u::Int )::Array
		
    visited[u] = true
    count += 1
    d[u] = count
    low[u] = d[u]
    for v in VV[u]
        if ! visited[v]
            push!(stack, [(u,v)])
            parent[v] = u
            DFV_visit( VV,out,count,visited,parent,d,low,stack, v )
            if low[v] >= d[u]
                push!(out, [outputComp(stack,u,v)])
            end
            low[u] = min( low[u], low[v] )
        else
            if ! (parent[u]==v) && (d[v] < d[u])
                push!(stack, [(u,v)])
            end
            low[u] = min( low[u], d[v] )
        end
    end
    out
end

# //////////////////////////////////////////////////////////////////////////////

function biconnectedComponent(EV)
    V = 1:maximum(EV)
    count = 0
    stack,out = [],[]
    visited = [false for v in V]
    parent = Union{Int, Array{Any,1}}[[] for v in V]
    d = Any[0 for v in V]
    low = Any[0 for v in V]    
    VV = Lar.verts2verts(EV)
    out = Any[]
    for u in V 
        if ! visited[u] 
            out = DFV_visit( VV,out,count,visited,parent,d,low,stack, u )
        end
    end
    out = [component for component in out if length(component) >= 1]
    EVs = [[map(sort∘collect,edges) for edges in union(comp...) if length(edges)>1] 
    		for comp in out]
    EVs = union(filter(x->!isempty(x), EVs)...)
    EVs = sort(EVs, lt=(x,y)->length(x)>length(y))
    return EVs
end

# //////////////////////////////////////////////////////////////////////////////

""" Algorithm for search of maximal biconnected components in a graph (no pao) """
function biconnected_components(EV::ChainOp)
    ps = Array{Tuple{Int, Int, Int}, 1}()
    es = Array{Tuple{Int, Int}, 1}()
    todel = Array{Int, 1}()
    visited = Array{Int, 1}()
    bicon_comps = Array{Array{Int, 1}, 1}()
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


#///// 9.2 Chain Complex Congruence  /////////////////////

#///// 9.2.1 Direct sum of local complexes  /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Merge congruent vertices and edges in `V` and `EV` """
function merge_vertices!(V::Points, EV::ChainOp, edge_map, err=1e-4)
    vertsnum = size(V, 1)
    edgenum = size(EV, 1)
    newverts = zeros(Int, vertsnum)
    # KDTree constructor needs an explicit array of Float64
    V = Array{Float64,2}(V)
    kdtree = KDTree(permutedims(V))

    # merge congruent vertices
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

""" Merge a local 2D complex; used in `planar_arrangement_2`, after componentgraph`` """
# codice incomprensibile, ma funzionante
function cell_merging(n, containment_graph, V, EVs, boundaries, shells, shell_bboxes)

   function bboxes(V::Points, indexes::ChainOp)
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

#///// 9.2.2 Mapping direct-sum to single complex   /////////////////////

#///// 9.3 Sparse matrix implementation  /////////////////////

#///// 9.3.1 The Graph-BLAS standard   /////////////////////

#///// 9.3.2 Using Julia Sparse Arrays   /////////////////////

#///// 9.3.3 Julia SuiteSparseGraphBLAS.jl   /////////////////////

#///// 10 Topological Gift Wrapping  /////////////////////

#///// 10.1 Meaning of chain operators  169  /////////////////////

#///// 10.1.1 Linear operators between linear spaces   /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Sparse binary matrix constructor for `Lar` d-cells, any `d`. """
function characteristicMatrix( CV::Cells )::ChainOp
	I,J,V = Int64[],Int64[],Int8[]
	for f=1:length(CV)
		for k in CV[f]
		push!(I,f)
		push!(J,k)
		push!(V,1)
		end
	end
	M = sparse(I,J,V)
	return M
end

# //////////////////////////////////////////////////////////////////////////////

""" Alternate method of function `u_coboundary_1` from `Cells` """
function u_coboundary_1( FV::Cells, EV::Cells, convex=true::Bool)::ChainOp
	copFV = characteristicMatrix(FV)
	copEV = characteristicMatrix(EV)
	out = u_coboundary_1( copFV::ChainOp, copEV::ChainOp, convex::Bool)
	return out
end

# //////////////////////////////////////////////////////////////////////////////

""" Unsigned function `coboundary_1` """
function u_coboundary_1( copFV::ChainOp, copEV::ChainOp, convex=true::Bool)::ChainOp
	temp = copFV * copEV'
	I,J,Val = Int64[], Int64[], Int8[]
	for j=1:size(temp,2)
		for i=1:size(temp,1)
			if temp[i,j] == 2
				push!(I,i)
				push!(J,j)
				push!(Val,1)
			end
		end
	end
	copFE = SparseArrays.sparse(I,J,Val)
	if !convex
		copFE = fix_redundancy(copFE,copFV,copEV)
	end
	return copFE
end


# //////////////////////////////////////////////////////////////////////////////

""" Alternate method of function `coboundary_1` from `Cells` """
function coboundary_1(FV::Cells, EV::Cells) 
	copFV = lar2cop(FV)
	I,J,Val = findnz(lar2cop(EV))
	copVE = sparse(J,I,Val)
	triples = hcat([[i,j,1]  for (i,j,v)  in zip(findnz(copFV * copVE)...) if v==2]...)
	I,J,Val = triples[1,:], triples[2,:], triples[3,:]
	Val = convert(Array{Int8,1},Val)
	copFE = sparse(I,J,Val)
	return copFE
end

# //////////////////////////////////////////////////////////////////////////////

""" Alternate method of function `coboundary_1` from `ChainOp` """
function coboundary_1( V::Points, FV::Cells, EV::Cells; convex=true::Bool, exterior=false::Bool)::ChainOp
	# generate unsigned operator's sparse matrix
	copFV = characteristicMatrix(FV)
	copEV = characteristicMatrix(EV)
   # greedy generation of incidence signs
   copFE = coboundary_1( V, copFV::ChainOp, copEV::ChainOp, convex, exterior )
	return copFE
end

# //////////////////////////////////////////////////////////////////////////////

""" Signed `coboundary_1` function, considering the LAR uncompleteness case """
function coboundary_1( V::Points, copFV::ChainOp, copEV::ChainOp, convex=true::Bool, exterior=false::Bool )::ChainOp

	copFE = u_coboundary_1( copFV::ChainOp, copEV::ChainOp, convex)
	EV = [findnz(copEV[k,:])[1] for k=1:size(copEV,1)]
	copEV = sparse(coboundary_0( EV::Cells ))
	for f=1:size(copFE,1)
		chain = findnz(copFE[f,:])[1]	#	dense
		cycle = spzeros(Int8,copFE.n)	#	sparse

		edge = findnz(copFE[f,:])[1][1]; sign = 1
		cycle[edge] = sign
		chain = setdiff( chain, edge )
		while chain != []
			boundary = sparse(cycle') * copEV
			_,vs,vals = findnz(dropzeros(boundary))

			rindex = vals[1]==1 ? vf = vs[1] : vf = vs[2]
			r_boundary = spzeros(Int8,copEV.n)	#	sparse
			r_boundary[rindex] = 1
			r_coboundary = copEV * r_boundary
			r_edge = intersect(findnz(r_coboundary)[1],chain)[1]
			r_coboundary = spzeros(Int8,copEV.m)	#	sparse
			r_coboundary[r_edge] = EV[r_edge][1]<EV[r_edge][2] ? 1 : -1

			lindex = vals[1]==-1 ? vi = vs[1] : vi = vs[2]
			l_boundary = spzeros(Int8,copEV.n)	#	sparse
			l_boundary[lindex] = -1
			l_coboundary = copEV * l_boundary
			l_edge = intersect(findnz(l_coboundary)[1],chain)[1]
			l_coboundary = spzeros(Int8,copEV.m)	#	sparse
			l_coboundary[l_edge] = EV[l_edge][1]<EV[l_edge][2] ? -1 : 1

			if r_coboundary != -l_coboundary  # false iff last edge
				# add edge to cycle from both sides
				rsign = rindex == EV[r_edge][1] ? 1 : -1
				lsign = lindex == EV[l_edge][2] ? -1 : 1
				cycle = cycle + rsign * r_coboundary + lsign * l_coboundary
			else
				# add last (odd) edge to cycle
				rsign = rindex==EV[r_edge][1] ? 1 : -1
				cycle = cycle + rsign * r_coboundary
			end
			chain = setdiff(chain, findnz(cycle)[1])
		end
		for e in findnz(copFE[f,:])[1]
			copFE[f,e] = cycle[e]
		end
	end
	if exterior && size(V,1)==2
		# put matrix in form: first row outer cell; with opposite sign )
		V = convert(Array{Float64,2},transpose(V))
		EV = convert(ChainOp, SparseArrays.transpose(boundary_1(EV)))

		outer = get_external_cycle(V::Points, copEV::ChainOp,
			copFE::ChainOp)
		copFE = [ -copFE[outer:outer,:];  copFE[1:outer-1,:];  copFE[outer+1:end,:] ]
		# induce coherent orientation of matrix rows (see examples/orient2d.jl)
		for k=1:size(copFE,2)
			spcolumn = findnz(copFE[:,k])
			if sum(spcolumn[2]) != 0
				row = spcolumn[1][2]
				sign = spcolumn[2][2]
				copFE[row,:] = -sign * copFE[row,:]
			end
		end
		return copFE
	else
		return copFE
	end
end

# //////////////////////////////////////////////////////////////////////////////


""" Utility function for `coboundary_0` """
function boundary_1( EV::Cells )::ChainOp
	out = characteristicMatrix(EV)'
	for e = 1:length(EV)
		out[EV[e][1],e] = -1
	end
	return out
end

""" `coboundary_0` (i.e. signed `EV`) function from `EV::Cells` """
coboundary_0(EV::Cells) = convert(ChainOp, transpose(boundary_1(EV::Cells)))


#///// 10.1.2 C3 space and δ2 :C2 → C3 matrix  /////////////////////

""" Main function of arrangement pipeline """
function space_arrangement(V::Points, EV::ChainOp, FE::ChainOp)
   fs_num = size(FE, 1)
   # strange but necessary cycle of computations to get FV::Cells algebraically
   FV = (abs.(FE) * abs.(EV)).÷ 2;
   FV = convert(ChainOp, FV);
   model = V', cop2lar(FV);
   sp_idx = spaceindex(model)

   rV = Points(undef, 0,3)
   rEV = SparseArrays.spzeros(Int8,0,0)
   rFE = SparseArrays.spzeros(Int8,0,0)

	depot_V = Array{Array{Float64,2},1}(undef,fs_num)
	depot_EV = Array{ChainOp,1}(undef,fs_num)
	depot_FE = Array{ChainOp,1}(undef,fs_num)
   for sigma in 1:fs_num
     print(sigma, "/", fs_num, "\r")
     nV, nEV, nFE = frag_face( V, EV, FE, sp_idx, sigma)
     depot_V[sigma] = nV
     depot_EV[sigma] = nEV
     depot_FE[sigma] = nFE
   end
	rV = vcat(depot_V...)
	rEV = SparseArrays.blockdiag(depot_EV...)
	rFE = SparseArrays.blockdiag(depot_FE...)

   rV, rcopEV, rcopFE = merge_vertices(rV, rEV, rFE)
   rcopCF = build_copFC(rV, rcopEV, rcopFE)
   return rV, rcopEV, rcopFE, rcopCF
end

#///// 10.2 Topological Gift Wrapping (TGW)   /////////////////////
#///// 10.2.1 Dim-independent TGW Algorithm   /////////////////////
#///// 10.2.2 2D Simplified version   /////////////////////


#///// 10.2.3 3D Standard version   /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" TGW algorithm implementation (pao) """
function build_copFC(rV, rcopEV, rcopFE)
#function build_copFC(V,FV,EV,copFE)

	# G&F -> Pao data structures
	V = convert(Points, rV')
	EV = cop2lar(rcopEV)
	fe = cop2lar(rcopFE)
	fv = [union([EV[e] for e in fe[f]]...) for f=1:length(fe)]
	FV = convert(Cells, fv)
	copFE = rcopFE    # alias useful for algorithm testing 
   VV = [[v] for v=1:size(V,2)]
   model = (V, [VV,EV,FV])
#Plasm.View(Plasm.numbering(.25)(model))

	copEF = transpose(copFE)
	FE = [SparseArrays.findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
	# Initializations
	m,n = size(copEF)
	#@show m,n;
	marks = zeros(Int,n);
	I = Int[]; J = Int[]; W = Int[];
	jcol = 0
	choose(marks) = findfirst(x -> x<2, marks)

	# Main loop (adding one copFC's column stepwise)
	while sum(marks) < 2n # no robust condition ... make better
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
		while nnz(cd2)≠0
			corolla = sparsevec([], Int[], m)
			# for each “hinge” τ cell
			for τ ∈ (.*)(SparseArrays.findnz(cd2)...)
				#compute the  coboundary
				#@show τ;
				tau = sparsevec([abs(τ)], Int[sign(τ)], m)  # ERROR: index out of bound here!
				bd1 = transpose(transpose(tau) * copEF)
				cells2D = SparseArrays.findnz(bd1)[1]
				# compute the  support
				inters = intersect(cells2D, SparseArrays.findnz(cd1)[1])
				if inters ≠ []
					pivot = inters[1]
					#@show pivot;
				else
				    #@show marks;
					error("no pivot")
				end
				# compute the new adj cell
				fan = ord(abs(τ),bd1,V,FV,EV,FE) # ord(pivot,bd1)
				if τ > 0
					adj = mynext(fan,pivot)
				elseif τ < 0
					adj = myprev(fan,pivot)
				end
				# orient adj
				if copEF[abs(τ),adj] ≠ copEF[abs(τ),pivot]
					corolla[adj] = cd1[pivot]
				else
					corolla[adj] = -(cd1[pivot])
				end
				#@show corolla;
			end
			# insert corolla cells in current cd1
			for (k,val) in zip(SparseArrays.findnz(corolla)...)
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
		append!(I,rows)
		append!(J,[ jcol for k=1:nnz(cd1) ])
		append!(W,vals)
	end
	copCF = sparse(J,I,W)
	return copCF
end



#///// 10.3 Chain space structure   /////////////////////

#///// 10.3.1 Chains, Cycles, and Boundaries   /////////////////////
#///// 10.3.2 Multi-cycle Boundaries Compute   /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Utility function for TGW in 3D """
function mynext(cycle, pivot)
	len = length(cycle)
	ind = findall(x -> x==pivot, cycle)[1]
	nextIndex = ind==len ? 1 : ind+1
	return cycle[nextIndex][1]
end


""" Utility function for TGW in 3D """
function myprev(cycle, pivot)
	len = length(cycle)
	ind = findall(x->x==pivot, cycle)[1]
	nextIndex = ind==1 ? len : ind-1
	return cycle[nextIndex][1]
end


#///// 10.3.3 Export B-rep to Graphics Files  /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Correct ordering of triangles (hence faces) about each boundary edge """
function ordering(triangles,V)
	normals = []
	v1,v2,v3 = triangles[1]
	if v1>v2 v1,v2 = v2,v1 end
	e3 = LinearAlgebra.normalize(V[:,v2]-V[:,v1])
	e1 = LinearAlgebra.normalize(V[:,v3]-V[:,v1])
	e2 = LinearAlgebra.normalize(cross(e1,e3))
	basis = [e1 e2 e3]
	transform = inv(basis)

	angles = []
	for (v1,v2,v3) in triangles
		w1 = LinearAlgebra.normalize(V[:,v3]-V[:,v1])
		w2 = transform * w1
		w3 = cross([0,0,1],w2)
		push!(normals,w3)
	end
	for k=1:length(normals)
		angle = atan(normals[k][2],normals[k][1])
		push!(angles,angle)
	end
	pairs = sort(collect(zip(angles,1:length(triangles))))
	order = [k for (angle,k) in pairs]
#@show order
	return order
end
# //////////////////////////////////////////////////////////////////////////////

""" From Lar topology to cells (1D chains, 2D chains, breps of 3D chains) """
function pols2tria(W, copEV, copFE, copCF) # W by columns
   V = convert(Points,W')
   #triangulated_faces = mytriangulate(V, [copEV, copFE])
   triangulated_faces = mytriangulate(V, [copEV, copFE])
   EVs = FV2EVs(copEV, copFE) # polygonal face fragments
   # triangulated_faces = [ item for (k,item) in enumerate(triangulated_faces)
   # 	if isdefined(triangulated_faces,k) && item ≠ Any[]  ]
   FVs = convert(Array{Cells}, triangulated_faces)
   CVs = []
   for cell in 1:copCF.m
      obj = []
          for f in copCF[cell, :].nzind
              triangles = triangulated_faces[f]
         append!(obj, triangles)
          end
      push!(CVs,obj)
      end
   V = convert(Points,V')
   return V,CVs,FVs,EVs
end
# //////////////////////////////////////////////////////////////////////////////

""" map 3D faces in z=0 and triangulate them; return a 3D Vector{Cells}  """
function mytriangulate(V::Points, cc::ChainComplex)   
   copEV, copFE = cc[1:2]

   triangulated_faces = Vector{Any}(undef, copFE.m)

   for f in 1:copFE.m
      if f % 10 == 0 print(".") end
      edges_idxs = copFE[f, :].nzind
      edge_num = length(edges_idxs)
      edges = zeros(Int, edge_num, 2)
      fv, edges = vcycle(copEV, copFE, f)
      if fv ≠ []
         vs = V[fv, :]
         v1 = LinearAlgebra.normalize(vs[2, :] - vs[1, :])
         v2 = [0, 0, 0]
         v3 = [0, 0, 0]
         err = 1e-8
         i = 3
         while -err < LinearAlgebra.norm(v3) < err
            v2 = LinearAlgebra.normalize(vs[i, :] - vs[1, :])
            v3 = LinearAlgebra.cross(v1, v2)
            i = i % size(vs,1) + 1
         end
         M = reshape([v1; v2; v3], 3, 3)
         vs = (vs*M)[:, 1:2]
         v = convert(Points, vs'[1:2,:])
         vmap = Dict(zip(fv,1:length(fv))) # vertex map
         mapv = Dict(zip(1:length(fv),fv)) # inverse vertex map
         trias = triangulate2d(v,edges)  # single face f
         triangulated_faces[f] = [[mapv[v] for v in tria] for tria in trias]
      end
   end
   return convert(Vector{Cells}, triangulated_faces)
end

# //////////////////////////////////////////////////////////////////////////////

""" Coherently orient the edges of f face """
function vcycle( copEV::ChainOp, copFE::ChainOp, f::Int )
  edges,signs = findnz(copFE[f,:])
  vpairs = [s>0 ? findnz(copEV[e,:])[1] :
              reverse(findnz(copEV[e,:])[1])
           for (e,s) in zip(edges,signs)]
  a = [pair for pair in vpairs if length(pair)==2]
  function mycat(a::Cells)
     out=[]
     for cell in a append!(out,cell) end
     return out
  end
  vs = collect(Set(mycat(a)))
  vdict = Dict(zip(vs,1:length(vs)))
  edges = [[vdict[pair[1]], vdict[pair[2]]] for pair in vpairs if length(pair)==2]
   return vs, edges
end

# //////////////////////////////////////////////////////////////////////////////

""" CDT Constrained Delaunay Triangulation """
function triangulate2d(V::Points, EV::Cells)
   points = convert(Array{Float64,2}, V')
	points_map = Array{Int,1}(collect(1:1:size(points)[1]))
   edges_list = convert(Array{Int,2}, hcat(EV...)')
   edge_boundary = [true for k=1:size(edges_list,1)] ## dead code !!
	#trias = constrained_triangulation2D(V::Lar.Points, EV::Lar.Cells)

 	trias = Triangle.constrained_triangulation(points,points_map,edges_list)
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


""" CDT Constrained Delaunay Triangulation """
function constrained_triangulation2D(V::Points, EV::Cells)
	triin = Triangulate.TriangulateIO()  # object generation
	triin.pointlist = V
	triin.segmentlist = hcat(EV...)
	(triout, vorout) = Triangulate.triangulate("pQ", triin) # exec triangulation
	trias = Vector{Int64}[c[:] for c in eachcol(triout.trianglelist)]
	return trias
end


#///// 11 Atoms & Generators of Solid Algebra   /////////////////////


#///// 11.1 Point-Polygon Classification   /////////////////////
# //////////////////////////////////////////////////////////////////////////////

""" Test of point containment in a polygon face """
function point_in_face(point, V::Points, copEV::ChainOp)

         function pointInPolygonClassification(V,EV)

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
           tiles = [[9,1,5],[8,0,4],[10,2,6]]
           b1,b2,b3,b4 = box
           """ code point position w.r.t query box using Bitwise OR """
           function tileCode(point)
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

      """ partial function; compute point classification w.r.t polygon edges """
        function pointInPolygonClassification0(pnt)
            x,y = pnt
            xmin,xmax,ymin,ymax = x,x,y,y
            tilecode = setTile([ymax,ymin,xmax,xmin])
            count,status = 0,0

            for k in 1:EV.m # loop on polygon edges
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
#    inner = !(pointInPolygonClassification(V, copEV)(point) == "p_out")
#    return inner
end

# //////////////////////////////////////////////////////////////////////////////

""" Test the whole arrangement pipeline; return the triangulated geometry """
function testarrangement(V,CV,FV,EV)
		cop_EV = coboundary_0(EV::Cells)
		cop_FE = convert(ChainOp, lar2cop(FV) * lar2cop(EV)'.÷ 2)
		W = convert(Points, V');
      #@show W, cop_EV, cop_FE;
		V, copEV, copFE, copCF = space_arrangement(
				W::Points, cop_EV::ChainOp, cop_FE::ChainOp);
      #@show V, copEV, copFE;
		V = convert(Points, V');
		V,CVs,FVs,EVs = pols2tria(V, copEV, copFE, copCF); ########### !!!
end;

# //////////////////////////////////////////////////////////////////////////////

V = [0.9101582 0.6379717 0.6592827 0.3870961 1.5062844 1.2340979 1.2554089 0.9832224 0.0843469 -0.0881351 0.2383533 0.0658713 0.3897303 0.2172483 0.5437367 0.3712547 -0.036538 -0.1151385 0.0390101 -0.0395905 0.2265768 0.1479763 0.3021248 0.2235243 0.7490819 0.6994406 0.6997791 0.6501378 1.1691151 1.1194738 1.1198123 1.070171 0.072093 0.3745675 -0.0689181 0.2335564 0.1464 0.4488745 0.0053889 0.3078634 0.1641349 0.1334626 0.1944844 0.1638121 0.3716781 0.3410058 0.4020276 0.3713553 0.0803545 0.6958805 0.147548 0.7630741 0.0877337 0.7032598 0.1549272 0.7704533 -0.1387411 -0.4444038 0.134104 -0.1715587 0.4015458 0.0958831 0.6743909 0.3687282 -0.1846788 -0.1108464 -0.2580307 -0.1841983 0.4551466 0.528979 0.3817947 0.4556271 0.304116 -0.1708791 0.0209969 -0.4539983 0.5142812 0.0392861 0.2311621 -0.2438331; 0.4028344 0.1519588 1.039914 0.7890385 0.5563971 0.3055215 1.1934767 0.9426012 0.022095 0.1761015 0.3624505 0.516457 -0.0625644 0.091442 0.277791 0.4317975 -0.0292294 0.0463186 0.2398719 0.3154199 -0.0839278 -0.0083798 0.1851735 0.2607215 0.6842707 0.6349679 1.1049786 1.0556758 0.7278259 0.6785231 1.1485337 1.0992309 0.5895819 0.4485708 0.4271526 0.2861415 0.8553432 0.7143321 0.6929139 0.5519028 0.7177194 0.748069 0.9259049 0.9562544 0.6917613 0.7221108 0.8999467 0.9302962 0.6810943 0.7482878 0.0802826 0.1474761 0.5471088 0.6143023 -0.0537029 0.0134906 0.5269135 0.7997587 1.1293121 1.4021573 0.3770618 0.6499069 0.9794604 1.2523055 0.6535834 0.5802315 1.292451 1.2190991 0.7352899 0.661938 1.3741575 1.3008056 0.2619841 -0.0211351 0.7783923 0.4952731 0.3177727 0.0346535 0.8341809 0.5510618; -0.1870511 0.4090752 -0.0334884 0.5626379 0.1497612 0.7458875 0.3033239 0.8994502 -0.0859219 0.2194615 -0.1705814 0.134802 0.1292543 0.4346377 0.0445948 0.3499782 0.4366695 0.6997843 0.3819711 0.6450859 0.5309756 0.7940904 0.4762772 0.739392 0.3690667 0.7890998 0.4126218 0.832655 0.4238204 0.8438536 0.4673756 0.8874087 0.6111982 0.6855052 0.8769595 0.9512665 0.8130543 0.8873613 1.0788156 1.1531226 0.3601404 0.5676836 0.3341822 0.5417255 0.3946086 0.6021518 0.3686504 0.5761937 0.1449814 0.1523606 0.0109959 0.0183751 0.7494938 0.756873 0.6155083 0.6228875 -0.07578 0.4645069 -0.2256318 0.3146551 0.3055579 0.8458448 0.1557061 0.695993 0.5562169 1.1960423 0.6379235 1.2777488 0.4917517 1.1315771 0.5734582 1.2132836 -0.2326003 -0.0224351 -0.1768117 0.0333535 0.3175492 0.5277144 0.3733378 0.583503]
CV = [[1, 2, 3, 4, 5, 6, 7, 8], [9, 10, 11, 12, 13, 14, 15, 16], [17, 18, 19, 20, 21, 22, 23, 24], [25, 26, 27, 28, 29, 30, 31, 32], [33, 34, 35, 36, 37, 38, 39, 40], [41, 42, 43, 44, 45, 46, 47, 48], [49, 50, 51, 52, 53, 54, 55, 56], [57, 58, 59, 60, 61, 62, 63, 64], [65, 66, 67, 68, 69, 70, 71, 72], [73, 74, 75, 76, 77, 78, 79, 80]]
FV = [[1, 2, 3, 4], [5, 6, 7, 8], [1, 2, 5, 6], [3, 4, 7, 8], [1, 3, 5, 7], [2, 4, 6, 8], [9, 10, 11, 12], [13, 14, 15, 16], [9, 10, 13, 14], [11, 12, 15, 16], [9, 11, 13, 15], [10, 12, 14, 16], [17, 18, 19, 20], [21, 22, 23, 24], [17, 18, 21, 22], [19, 20, 23, 24], [17, 19, 21, 23], [18, 20, 22, 24], [25, 26, 27, 28], [29, 30, 31, 32], [25, 26, 29, 30], [27, 28, 31, 32], [25, 27, 29, 31], [26, 28, 30, 32], [33, 34, 35, 36], [37, 38, 39, 40], [33, 34, 37, 38], [35, 36, 39, 40], [33, 35, 37, 39], [34, 36, 38, 40], [41, 42, 43, 44], [45, 46, 47, 48], [41, 42, 45, 46], [43, 44, 47, 48], [41, 43, 45, 47], [42, 44, 46, 48], [49, 50, 51, 52], [53, 54, 55, 56], [49, 50, 53, 54], [51, 52, 55, 56], [49, 51, 53, 55], [50, 52, 54, 56], [57, 58, 59, 60], [61, 62, 63, 64], [57, 58, 61, 62], [59, 60, 63, 64], [57, 59, 61, 63], [58, 60, 62, 64], [65, 66, 67, 68], [69, 70, 71, 72], [65, 66, 69, 70], [67, 68, 71, 72], [65, 67, 69, 71], [66, 68, 70, 72], [73, 74, 75, 76], [77, 78, 79, 80], [73, 74, 77, 78], [75, 76, 79, 80], [73, 75, 77, 79], [74, 76, 78, 80]]
EV = [[1, 2], [3, 4], [5, 6], [7, 8], [1, 3], [2, 4], [5, 7], [6, 8], [1, 5], [2, 6], [3, 7], [4, 8], [9, 10], [11, 12], [13, 14], [15, 16], [9, 11], [10, 12], [13, 15], [14, 16], [9, 13], [10, 14], [11, 15], [12, 16], [17, 18], [19, 20], [21, 22], [23, 24], [17, 19], [18, 20], [21, 23], [22, 24], [17, 21], [18, 22], [19, 23], [20, 24], [25, 26], [27, 28], [29, 30], [31, 32], [25, 27], [26, 28], [29, 31], [30, 32], [25, 29], [26, 30], [27, 31], [28, 32], [33, 34], [35, 36], [37, 38], [39, 40], [33, 35], [34, 36], [37, 39], [38, 40], [33, 37], [34, 38], [35, 39], [36, 40], [41, 42], [43, 44], [45, 46], [47, 48], [41, 43], [42, 44], [45, 47], [46, 48], [41, 45], [42, 46], [43, 47], [44, 48], [49, 50], [51, 52], [53, 54], [55, 56], [49, 51], [50, 52], [53, 55], [54, 56], [49, 53], [50, 54], [51, 55], [52, 56], [57, 58], [59, 60], [61, 62], [63, 64], [57, 59], [58, 60], [61, 63], [62, 64], [57, 61], [58, 62], [59, 63], [60, 64], [65, 66], [67, 68], [69, 70], [71, 72], [65, 67], [66, 68], [69, 71], [70, 72], [65, 69], [66, 70], [67, 71], [68, 72], [73, 74], [75, 76], [77, 78], [79, 80], [73, 75], [74, 76], [77, 79], [78, 80], [73, 77], [74, 78], [75, 79], [76, 80]]

GL.VIEW([ GL.GLPol(V,CV, GL.COLORS[1], 0.5) ]);

V,CVs,FVs,EVs = testarrangement(V,CV,FV,EV);

GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2,99,1));
GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
sort!( CVs; by=x->length(x) );
GL.VIEW(GL.GLExplode(V,CVs[1:end-2],4.,4.,4.,99,0.5));
GL.VIEW(GL.GLExplode(V,[CVs[end-1]],1.,1.,1.,99,0.5));


# //////////////////////////////////////////////////////////////////////////////

using Plasm
Cubes = STRUCT([ CUBE(1), T([1,2,3])([.5,.5,.5]), CUBE(1)])
Brep = ToBoundaryForm(Cubes)
V, FV = UKPOL(Brep)
W = hcat(V...)

function simplexfacets(simplices)
   @assert hcat(simplices...) isa Matrix
   out = Array{Int64,1}[]
	for simplex in simplices
		for v in simplex
			facet = setdiff(simplex,v)
			push!(out, facet)
		end
	end
	# remove duplicate faces
	return sort(collect(Set(out)))
end

EV = simplexfacets(FV)
_,CV,_ = UKPOL(Cubes)
V,CVs,FVs,EVs = testarrangement(W,CV,FV,EV);


#///// 11.2 Set-Merbership Classification (SMC)  /////////////////////
#///// 11.2.1 Point to Boolean Atom Mapping   /////////////////////
#///// 11.2.2 Structure of Variables as Binary Tuples   /////////////////////
#///// 11.3 Evaluation of CSG Expression  173 11.3.1 JuliaNativeBitwiseOperators174 11.3.2 Convert AST to set algebra expression  176 11.3.3 Bitwise resolution of set algebra expressions           /////////////////////
#///// 11.4 CSG to Boundary Compute 176 11.4.1 Apply a pipe of Chain Maps  178 11.4.2 Export to Graphics Formats   /////////////////////



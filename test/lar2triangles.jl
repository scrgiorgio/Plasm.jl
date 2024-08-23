
# //////////////////////////////////////////////////////////////////////////////
using Plasm


# //////////////////////////////////////////////////////////////////////////////

function constrained_triangulation2D(V::Points, EV::Cells)
	triin = Triangulate.TriangulateIO()    # object generation
	triin.pointlist = permutedims(V)
	triin.segmentlist = hcat(EV...)
	(triout, vorout) = Triangulate.triangulate("pQ", triin)  # exec triangulation
	trias = Array{Int64,1}[c[:] for c in eachcol(triout.trianglelist)]
	return trias
end

# //////////////////////////////////////////////////////////////////////////////
""" CDT Constrained Delaunay Triangulation """
function triangulate2d(V, EV) 
   points = convert(Matrix{Float64}, V')
#	points_map = Array{Int,1}(collect(1:1:size(points,2)[1]))
#   edges_list = convert(Array{Int,2}, hcat(EV...)')
#   edge_boundary = [true for k=1:size(edges_list,1)] ## dead code !!
#  trias = Triangle.constrained_triangulation(points,points_map,edges_list)
 	
	trias = constrained_triangulation2D(points, EV)
	innertriangles = Array{Int,1}[]
	copEV = lar2cop(EV)
	for (u,v,w) in trias
		point = (points[u,:]+points[v,:]+points[w,:])./3
		inner = point_in_face(point, points::Points, copEV::ChainOp)
		if inner
			push!(innertriangles,[u,v,w])
		end
	end
   return innertriangles
end

# //////////////////////////////////////////////////////////////////////////////
""" return ordered vertices  and edges of the 1-cycle f """
function cycle( EV, FE, f::Int )
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
   edges = [[a,b] for (a,b) in zip(ordered[1:end-1],ordered[2:end])]
   return Array{Int}(ordered[1:end-1]), edges
end

# //////////////////////////////////////////////////////////////////////////////
function lar2triangles(V, EV, FV, FE) 
   V = size(V,1)==3 ? permutedims(V) : V
   triangulated_faces = Vector{Any}(undef, length(FE))

   for f in 1:length(FE)
      edges_idxs = FE[f]
      edge_num = length(edges_idxs)
      fv, edges = cycle(EV, FE, f) 
      # look for independent vector triple
      points = V[fv,:]
      vmap = Dict(zip(fv,1:length(fv))) # inverse vertex map
      mapv = Dict(zip(1:length(fv),fv)) # inverse vertex map
      edges = [[vmap[A],vmap[B]] for (A,B) in edges]
      v1 = LinearAlgebra.normalize(points[2,:] - points[1,:])
      v2 = [0, 0, 0];   v3 = [0, 0, 0]
      err = 1e-8; i = 3 
      while -err < LinearAlgebra.norm(v3) < err
         v2 = LinearAlgebra.normalize(points[i,:] - points[1,:])
         v3 = LinearAlgebra.cross(v1, v2)
         i = i % size(points,1) + 1
      end
      # independent vector triple in face f 
      M = [v1 v2 v3] 
      projected = (points * M)[:,1:2]
      trias = triangulate2d(permutedims(projected),edges)  # single face f
      triangulated_faces[f] = [[mapv[v] for v in tria] for tria in trias]

   end
   return triangulated_faces
end
   
# ////// TEST /////////////////////////////////////////////////////////////   

using Plasm, SparseArrays,LinearAlgebra,Triangulate

function test_triangulation()
# dati all'uscita di "space_arrangement" dell'esempio atomsdrawing.jl
   (V, copEV, copFE, copCF) = ([0.70710678118655 0.70710678118655 0.0; 2.829225697485796e-17 0.0 0.0; -0.70710678118655 0.70710678118655 0.0; 0.0 1.4142135623731 0.0; 0.4142135623731001 1.0 0.0; 1.0683217716804808e-17 1.0 0.0; 0.70710678118655 0.70710678118655 1.0; 0.0 0.0 1.0; 0.0 1.4142135623731 1.0; 0.4142135623731001 1.0 1.0; -0.70710678118655 0.70710678118655 1.0; 0.0 1.0 1.0; 1.0 0.0 0.0; 1.0 1.0 0.0; 1.0 0.0 1.0; 1.0 1.0 1.0], sparse([1, 5, 10, 1, 2, 6, 9, 20, 2, 3, 16, 3, 4, 13, 4, 5, 7, 14, 21, 6, 7, 27, 8, 10, 11, 8, 9, 15, 18, 23, 12, 13, 17, 11, 12, 14, 19, 28, 15, 16, 17, 18, 19, 27, 20, 22, 24, 21, 22, 26, 23, 24, 25, 25, 26, 28], [1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16], Int8[-1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1], 28, 16), sparse([2, 3, 10, 1, 6, 1, 7, 1, 4, 2, 5, 10, 1, 2, 13, 1, 2, 15, 3, 8, 16, 3, 6, 11, 13, 3, 5, 5, 8, 16, 4, 9, 4, 7, 4, 5, 14, 15, 6, 9, 6, 7, 7, 9, 8, 9, 13, 8, 9, 15, 10, 11, 10, 14, 10, 12, 11, 16, 11, 12, 12, 16, 12, 14, 13, 15, 14, 16], [1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 10, 10, 11, 11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28], Int8[-1, 1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, -1, 1, 1, -1], 16, 28), sparse([2, 4, 2, 3, 1, 3, 2, 4, 1, 3, 2, 4, 2, 4, 2, 3, 2, 4, 1, 2, 1, 2, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2], [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16], [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1], 4, 16));
   
   EV = AA(sort)([findnz(copEV[k,:])[1] for k=1:copEV.m]); # vertices per edge
   FE = AA(sort)([findnz(copFE[k,:])[1] for k=1:copFE.m]); # edges per face
   FV = [union(CAT([EV[e]  for e in FE[f]])) for f=1:length(FE)]; # verts x face
   
   triangulated_faces = lar2triangles(V, EV, FV, FE)
end   
   
triangulated_faces = test_triangulation()
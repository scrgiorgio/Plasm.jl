
using Triangulate
export triangulate2d,find_cycle, lar2triangles

# //////////////////////////////////////////////////////////////////////////////
""" return ordered vertices  and edges of the 1-cycle f """
function find_cycle( EV, FE, f::Int )
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
""" input old LAR consistent data; output triangulated_faces """
function lar2triangles(V, EV, FV, FE) 
   V = size(V,1)==3 ? permutedims(V) : V
   triangulated_faces = Vector{Any}(undef, length(FE))

   for f in 1:length(FE)
      edges_idxs = FE[f]
      edge_num = length(edges_idxs)
      fv, edges = find_cycle(EV, FE, f) 
      # look for independent vector triple
      points = V[fv,:]
      vmap = Dict(zip(fv,1:length(fv))) # vertex map
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
   
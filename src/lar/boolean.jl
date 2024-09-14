
using Plasm
using Random
using SparseArrays
using LinearAlgebra

# //////////////////////////////////////////////////////////////////////////////
function random(a::Float64,b::Float64)
  return a+rand()*(b-a)
end

# //////////////////////////////////////////////////////////////////////////////
function random_dir()
  return normalized([random(-1.0,1.0) for I in 1:3])
end





# //////////////////////////////////////////////////////////////////////////////
function ray_face_intersection(ray_origin::Vector{Float64}, ray_dir::Vector{Float64}, lar::Lar, F::Int)

  face          = lar.C[:FV][F]
  face_points3d = lar.V[:,face]
	plane=plane_create(face_points3d)

	hit3d, t=plane_ray_intersection(ray_origin, ray_dir, plane)
	if isnothing(hit3d)
		return nothing,nothing
	end

	# i need to check if the hit is really inside the face
	#   to do so I project all in 2d and use the 2d classify point
	project = project_points3d(face_points3d; double_check=true) # scrgiorgio: remove double check 

	hit2d         = project(hit3d[:,:])[:,1]
	face_points2d = project(face_points3d)

  vdict=Dict(vertex_index => I for (I, vertex_index) in enumerate(face))
  face_edges=[ [vdict[a],vdict[b]] for (a,b) in lar.C[:EV] if a in face && b in face]
  # @show(face_edges)

	classify = classify_point(hit2d, BYROW(face_points2d), face_edges) 
  if classify == "p_out"
    return nothing,nothing
  end
  
  # println("ray ", ray_origin, " ",ray_dir, " intersecting  with face ",F, " ", BYROW(face_points3d),"---",hit2d, BYROW(face_points2d), face_edges)
  return hit3d, t

end
export ray_face_intersection


# //////////////////////////////////////////////////////////////////////////////
function is_internal_point(lar::Lar, ray_origin::Vector{Float64}, ray_dir::Vector{Float64})
  num_intersections=0
  for F in 1:length(lar.C[:FV])
    hit,t=ray_face_intersection(ray_origin, ray_dir, lar, F)
    if !isnothing(hit)
      num_intersections+=1
    end
  end
  return (num_intersections % 2)==1
end
export is_internal_point


# //////////////////////////////////////////////////////////////////////////////
function find_internal_point(lar::Lar;max_attempts=10) 

  # @show(lar)

  #VIEWCOMPLEX(lar)

	# TODO: reintroduce space index to avoid O(N^2) complexity

  b1,b2=lar_bounding_box(lar; only_used_vertices=true)
  
  # i want to move out of the complex
  move_out = LinearAlgebra.norm(b2 - b1)   

  #@show(b1,b2)
  #@show(move_out)

  for attempt in 1:max_attempts

    # choose a random point inside the box and a random direction
    inside_bbox = [random(b1[I],b2[I]) for I in 1:3 ]

    external_point = inside_bbox + move_out * random_dir()
    ray_dir    = normalized(inside_bbox-external_point)

    #@show("Trying", external_point, ray_dir)

    distances=Vector{Float64}()
    for F in 1:length(lar.C[:FV])
      hit,distance=ray_face_intersection(external_point, ray_dir, lar, F)
      if !isnothing(hit)
        push!(distances,distance)
        #@show("hit", external_point,ray_dir,distance)
      end
    end
    
    # I should be outside after I've hit all faces (starting from the outside)
    num_hits=length(distances)
    if (num_hits < 2) || (num_hits % 2) != 0
      continue
    end

    # I am going from OUT to in the OUT then in... finally OUT
    # I want to find the inner IN range which is bigger
    begin
      distances=sort(distances)
      #@show(distances)
      best_delta,internal_point=nothing,nothing
      for I in 1:length(distances)-1
        delta=distances[I+1]-distances[I]
        if ((I % 2) ==1) && delta>LAR_DEFAULT_ERR # only if is inside and delta is not too small
          distance=distances[I]+delta/2
          if isnothing(best_delta) || delta>best_delta
            internal_point=external_point+ray_dir*distance
            #@show("new best", external_point, ray_dir, distance,internal_point, best_delta,delta)
            best_delta=delta
          else
            #@show("not the best",delta)
          end
        end
      end
    end

    # try another time...
    if isnothing(internal_point)
      continue
    end

    # how to go from internal_point to outside
    ray_dir=-ray_dir 
    internal_external_distance=LinearAlgebra.norm(external_point-internal_point)
    distances=[(internal_external_distance-distance) for distance in distances]
    @assert(is_internal_point(lar, internal_point, ray_dir))
    #@show(internal_point, ray_dir, distances)
    return internal_point, ray_dir, distances

  end

	error("cannot find internal point")
end
export find_internal_point


function Union(v::AbstractArray)        return any(v) end
function Intersection(v::AbstractArray) return all(v) end
function Difference(v::AbstractArray)   return v[1] && !any(v[2:end]) end
function Xor(v::AbstractArray)          return (length([it for it in v if it]) % 2)==1  end

export Union, Intersection,Difference,Xor

# ////////////////////////////////////////////////////////////////
function get_outer_atom(atoms)
  diags =[LinearAlgebra.norm(b[2] - b[1]) for b in [lar_bounding_box(atom, only_used_vertices=true) for atom in atoms]]
  outer, outer_index = findmax(diags)
  return outer, outer_index
end
export get_outer_atom

# ////////////////////////////////////////////////////////////////
function bool3d(arrangement::Lar; input_args=[], bool_op=Union, debug_mode=true)

  atom_faces=arrangement.C[:CF]

  atoms=[]
  for sel in atom_faces
    atom=SELECT(arrangement, sel)
    push!(atoms,atom)
    # VIEWCOMPLEX(atom, explode=[1.4,1.4,1.4])
  end
  
  outer_atom, outer_index = get_outer_atom(atoms)
  deleteat!(atoms,      outer_index)
  deleteat!(atom_faces, outer_index)
  
  internal_points=[find_internal_point(atom) for atom in atoms] 
  
  bool_matrix=zeros(Bool,length(atoms),length(input_args))
  
  for (A,(atom, (internal_point, ray_dir, distances))) in enumerate(zip(atoms, internal_points))
    # @show(internal_point, ray_dir, distances)
  
    # find for each input arg if the internal point is inside or outside
    for (I,input_arg) in enumerate(input_args)
      bool_matrix[A,I]=is_internal_point(input_arg, internal_point, ray_dir)
    end
  
    if debug_mode
      points = GLBatch(POINTS); points.point_size = 8
      lines  = GLBatch(LINES) ; lines.line_width  = 2
      append!(points.colors.vector,RED); append!(points.vertices.vector, internal_point)
      for distance in distances
        append!(points.colors.vector, GREEN ); append!(points.vertices.vector, internal_point+distance*ray_dir)
        append!(lines.colors.vector,  YELLOW); append!(lines.vertices.vector, internal_point)
        append!(lines.colors.vector,  YELLOW); append!(lines.vertices.vector, internal_point+distance*ray_dir)
      end
  
      text=GLText(join([string(it) for it in bool_matrix[A,:]], " "), center=internal_point, color=BLACK)
      VIEWCOMPLEX(atom, show=["V", "EV"], explode=[1.0,1.0,1.0], user_batches=[points,lines,text...])
    end
  end
  
  sel=Cell()
  for (A,row) in enumerate(eachrow(bool_matrix))
    if bool_op(row)
      append!(sel,atom_faces[A])
    end
  end

  return SELECT(arrangement, remove_duplicates(sel))

end
export bool3d


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


# //////////////////////////////////////////////////////////////////////////////
function bool3d(assembly::Hpc, arrangement::Lar, CF::Cells; debug_mode=true)

  atoms=[]
  for sel in CF
    atom=SELECT(arrangement,sel)
    push!(atoms,atom)
  end

  # remove outer atoms (i.e. the one that contains all other atoms) by using diagonal measurement
  begin
    diags =[LinearAlgebra.norm(b[2] - b[1]) for b in [lar_bounding_box(atom,only_used_vertices=true) for atom in atoms]]
    __outer, outer_index = findmax(diags)
    deleteat!(atoms, outer_index)
  end

  # create the boolean matrix
  # TODO: probably if I work on the overall assembly I have to test few faces (ref. old Alberto's face span code)
  
  atoms=[ [atom,find_internal_point(atom)...]  for atom in atoms]

  # view atoms 
  if debug_mode
    for (atom, ray_origin, ray_dir) in atoms

      batches=BATCHES(atom)

      begin
        points = GLBatch(POINTS)
        points.point_size = 3
        points.point_color = RED
        append(points.vertices.vector, ray_origin)
        push!(batches,points)
      end

      begin
        lines = GLBatch(LINES)
        lines.line_width = 3
        lines.point_color = YELLOW
        append(lines.vertices.vector, (ray_origin+0.0*ray_dir))
        append(lines.vertices.vector, (ray_origin+1.0*ray_dir))
        push!(batches,points)
      end

      VIEWBATCHES(batches)

    end
  end

  # TOPOS is needed to isolate input arguments (will create one lar for each input arg)
  ret=[]
  args = [LAR(it) for it in TOPOS(assembly)]
  for (atom, ray_origin, ray_dir) in atoms
    push!(ret,[is_internal_point(arg, ray_origin, ray_dir) for arg in args])
  end

  return ret
end
export bool3d


using Plasm
using Random
using SparseArrays
using LinearAlgebra

# //////////////////////////////////////////////////////////////////////////////
function random(a::Float64,b::Float64)
  return a+rand()*(b-a)
end

function random_dir()
  return normalized([random(-1.0,1.0) for I in 1:3])
end


# //////////////////////////////////////////////////////////////////////////////
function bbox(lar::Lar; only_used_vertices=false)
  V=only_used_vertices ? lar.V[:, CAT(lar.C[:EV]) ] : lar.V
  return collect([vec(it) for it in bbox_create(V)])
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
  #@show(hit2d)
  #@show(BYROW(face_points2d))
  #@show(face_edges)
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
function find_internal_point(lar::Lar;num_attempts=13) 

  println("find_internal_point")
  @show(lar)

  #delete!(lar.text, :EV)
  #delete!(lar.text, :FV)
  # batches=LAR_BATCHES(lar)
  #View(batches)

	# TODO: reintroduce space index to avoid O(N^2) complexity

  b1,b2=bbox(lar; only_used_vertices=true)
  
  # i want to move out of the complex
  move_out = LinearAlgebra.norm(b2 - b1)   

  @show(b1,b2)
  @show(move_out)

  for attempt in 1:num_attempts

    # choose a random point inside the box and a random direction
    inside_bbox = [random(b1[I],b2[I]) for I in 1:3 ]

    external_point = inside_bbox + move_out * random_dir()
    ray_dir    = normalized(inside_bbox-external_point)

    @show("Trying", external_point, ray_dir)

    ts=[]
    for F in 1:length(lar.C[:FV])
      hit,t=ray_face_intersection(external_point, ray_dir, lar, F)
      if !isnothing(hit)
        push!(ts,t)
        @show("hit", external_point,ray_dir,t)
      end
    end
    
    # I should be outside after I've hit all faces (starting from the outside)
    num_hits=length(ts)
    if (num_hits < 2) || (num_hits % 2) != 0
      continue
    end

    # take the internal range which is bigger
    ts=sort(ts)
    @show(ts)

    best_delta=nothing
    internal_point=nothing
    for I in 1:length(ts)-1
      delta=ts[I+1]-ts[I]
      if ((I % 2) ==1) && delta>LAR_DEFAULT_ERR # only if is inside and delta is not too small
        t=ts[I]+delta/2
        if isnothing(best_delta) || delta>best_delta
          internal_point=external_point+ray_dir*t
          @show("new best", external_point, ray_dir, t,internal_point, best_delta,delta)
          best_delta=delta
        else
          @show("not the best",delta)
        end
      end
    end

    if isnothing(internal_point)
      return nothing
    end

    # how to go from internal_point to outside
    ray_dir=-ray_dir
    @assert(is_internal_point(lar, internal_point, ray_dir))
    @show(internal_point, ray_dir)
    return internal_point, ray_dir 
  end

  aaa()
	error("cannot find internal point")
end
export find_internal_point

# ///////////////////////////////////////////////////////////
"""
function example_save_points()

  ray_dir=random_dir()

  @show(is_internal_point(atoms[1], [1.0,1.0,-2.0], ray_dir))
  aaa()

  points=[]
  for z in range(-2.0,+2.0,step=0.05)
    for y in range(-2.0,+2.0,step=0.05)
      for x in range(-2.0,+2.0,step=0.05)
        if is_internal_point(atoms[1], [x,y,z], [0.0,0.0,1.0])
          push!(points,[x,y,z])
        end
      end
    end
  end

  open("tmp.ply", "w") do file
    println(file, "ply")
    println(file, "format ascii 1.0")
    println(file, "element vertex ", length(points))
    println(file, "property float x")
    println(file, "property float y")
    println(file, "property float z")
    println(file, "property float intensity")
    println(file, "end_header")
    for (x,y,z) in points
      println(file,x," ",y," ",z," ", 1)
    end
  end
end
"""

# //////////////////////////////////////////////////////////////////////////////
function bool3d(assembly::Hpc, arrangement::Lar, CF::Cells; debug_mode=true)

  #@show(arrangement)
  #@show(CF)

  atoms=[]
  for sel in CF
    atom=SELECT_FACES(arrangement,sel)
    push!(atoms,atom)
  end

  # remove outer atoms (i.e. the one that contains all other atoms) by using diagonal measurement
  begin
    diags =[LinearAlgebra.norm(b[2] - b[1]) for b in [bbox(atom,only_used_vertices=true) for atom in atoms]]
    __outer, outer_index = findmax(diags)
    deleteat!(atoms, outer_index)
  end

  # todo: remove
  # atoms=[atoms[1]]


  # create the boolean matrix
  # TODO: probably if I work on the overall assembly I have to test few faces (ref. old Alberto's face span code)
  
  atoms=[[atom,find_internal_point(atom)...]  for atom in atoms]

  @show(atoms)

  # view atoms 
  if debug_mode
    for (atom,ray_origin,ray_dir) in atoms
      batches=LAR_BATCHES(atom)

      begin
        batch = GLBatch(POINTS)
        push!(batches,batch)
        batch.point_size = 4
        batch.point_color = RED
        push!(batch.vertices.vector, ray_origin...)
      end

      begin
        batch = GLBatch(LINES)
        push!(batches,batch)
        batch.line_width = 3
        batch.point_color = YELLOW
        push!(batch.vertices.vector, (ray_origin+0.0*ray_dir)...)
        push!(batch.vertices.vector, (ray_origin+1.0*ray_dir)...)
      end

      View(batches)
    end
  end

  ret=[]
  args = [LAR(it) for it in TOPOS(assembly)]
  for (atom, ray_origin, ray_dir) in atoms
    push!(ret,[is_internal_point(arg, ray_origin, ray_dir) for arg in args])
  end

  @show(ret)
  return ret
end
export bool3d

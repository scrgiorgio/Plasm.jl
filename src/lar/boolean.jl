
using Plasm
using Random
using SparseArrays
using LinearAlgebra

# //////////////////////////////////////////////////////////////////////////////
function random_float(a::Float64,b::Float64)
  return a+rand()*(b-a)
end
export random_float

# //////////////////////////////////////////////////////////////////////////////
function random_dir(dim::Int)
  return normalized([random_float(-1.0,1.0) for I in 1:dim])
end
export random_dir

# ////////////////////////////////////////////////////////////////////////
function lar_connected_components(seeds::Cell, get_connected::Function)::Cells
  ret = []
  assigned = Set()

  function visit(component, cur::Int)
    if cur in assigned return end
    push!(assigned, cur)
    @assert !(cur in component)
    push!(component, cur)
    for other in get_connected(cur)
      visit(component, other)
    end
  end

  for seed in seeds
      if seed in assigned continue end
      component = Set()
      visit(component, seed)
      push!(ret, collect(component))
  end

  return ret
end

# //////////////////////////////////////////////////////////////////////////////
function guess_boundary_faces(lar::Lar, faces::Vector; max_attempts=1000)::Cell

  pdim=size(lar.V, 1)

  for attempt in 1:max_attempts
    b1,b2=lar_bounding_box(lar; only_used_vertices=true)
    move_out = pdim*LinearAlgebra.norm(b2 - b1)   
    inside_bbox = [random_float(b1[I],b2[I]) for I in 1:pdim ]
    external_point = inside_bbox + move_out * random_dir(pdim)
    ray_dir = normalized(inside_bbox-external_point)
    distances=[]
    for F in faces
      hit, distance=ray_face_intersection(external_point, ray_dir, lar, F)
      if !isnothing(hit)
        push!(distances,[distance,F])
      end
    end

    # this means I need two hit and should start outside and end in outside
    if length(distances) >=2 && (length(distances) % 2) == 0
      println("guess_boundary_faces #attempt=",attempt) # , " distances=", distances)
      distances=sort(distances)
      return [ distances[1][end], distances[end][end]]
    end

    # println("FAILED guess_boundary_faces #attempt=", attempt, " distances=", distances)

  end
  @assert(false)
end


# ////////////////////////////////////////////////////////////////
function SPLIT(lar::Lar; debug_mode=false)::Tuple{Lar,Lar}

  pdim=size(lar.V, 1)
  @assert(pdim==3)

  function compute_atom_bbox(lar::Lar, atom::Cell)
    points=PointsNd()
    for F in atom
      append!(points,[p for p in eachcol(lar.V[:, lar.C[:FV][F] ])])
    end
    return bbox_create(points)
  end

  atoms=[cf for cf in lar.C[:CF]]

  # connected atoms
  begin
    #for (A,atom) in enumerate(atoms)
    #  println("Atom ",A," ",atom)
    #end
    components=lar_connected_components(collect(eachindex(atoms)), A -> [B for B in eachindex(atoms) if A!=B && length(intersect(Set(atoms[A]),Set(atoms[B])))>0 ])
    components=[ [atoms[jt] for jt in it] for it in components]

    for (C,component) in enumerate(components)
      # not 100% sure about it
      if length(component)>2
        components[C]=remove_duplicates(component)
      else
        components[C]=[normalize_cell(it) for it in component]
      end
    end

    #for (C,component) in enumerate(components)
    #  println("Component ",C," ",component)
    #end
  end

  lar_outers=lar_copy(lar);lar_outers.C[:CF]=[]
  lar_inners=lar_copy(lar);lar_inners.C[:CF]=[]
  
  for (C,component) in enumerate(components)

    faces=remove_duplicates(vcat(component...))
    num_atoms=length(component)
    #println("# components ",C, " num_atoms=", num_atoms, " faces=",faces)
    #println(component)

    # there should be at least the internal and the external
    @assert(length(component)>=2) 

    outer, inners=nothing,[]

    # there is one outer cell, and one inner cell and they must be the same
    if num_atoms==2
      @assert(component[1]==component[2])
      println("length is 2, one outer, one inner")
      outer,inners=component[1],[component[2]]

    else
      
      # try guessing with the bounding box
      begin
        bboxes=[compute_atom_bbox(lar, atom) for (A,atom) in enumerate(component)]
        for (O,maybe_outer) in enumerate(component)
          is_outer=true
          maybe_inners=[]
          for (I,maybe_inner) in enumerate(component)
            if I!=O 
              push!(maybe_inners, maybe_inner)
              is_outer=is_outer && bbox_contain(bboxes[O], bboxes[I]) && bboxes[O]!=bboxes[I]
            end
          end
          if is_outer
            @assert(isnothing(outer))
            outer,inners=maybe_outer,maybe_inners
            println("Found outer by using bounding box")
          end
        end
      end

      # try guessing with boundary faces (can be slow and error-prone)
      if isnothing(outer)
        filtered=component
        while length(filtered)>1
          boundary_faces=guess_boundary_faces(lar, faces)
          filtered = [atom for atom in filtered if all([F in atom for F in boundary_faces])]
        end
        @assert(length(filtered)==1)
        outer =filtered[1]
        inners=[atom for atom in component if atom!=outer]
      end

    end

    # println("outer is ",outer)

    # topology check: 
    begin
      @assert(!isnothing(outer))
      @assert(length(inners)>=1)
      @assert((length(inners)+1)==length(component))

      # all outer faces should in any of the inners
      for F in outer 
        @assert(F in vcat(inners...))
      end
    end

    push!(lar_outers.C[:CF], outer)
    append!(lar_inners.C[:CF], inners)

  end

  return lar_outers, lar_inners

end

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

# //////////////////////////////////////////////////////////////////////////////
function ray_face_intersection(ray_origin::PointNd, ray_dir::PointNd, lar::Lar, idx::Int)

  pdim=size(lar.V, 1)

  if pdim===2

    # scrgiorgio: not tested. 
    a,b = lar.C[:EV][idx]
    P1  = lar.V[:, a]
    P2  = lar.V[:, b]

    v = P1 - ray_origin  
    w = P2 - P1 
    A = [ray_dir -w]
    if abs(det(A)) < LAR_DEFAULT_ERR 
      return nothing, nothing # parallel or coincidents
    end
    t_u = A \ v
    t, u = t_u[1], t_u[2]
    if t >= 0.0 && 0.0 <= u <= 1.0
        return ray_origin + t * ray_dir, t
    else
      # println("here ", t," ", u)
      return nothing,nothing
    end

  else
    @assert(pdim==3)

    face          = lar.C[:FV][idx]
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
    
    # println("ray ", ray_origin, " ",ray_dir, " intersecting  with face ",idx, " ", BYROW(face_points3d),"---",hit2d, BYROW(face_points2d), face_edges)
    return hit3d, t

  end

end

# //////////////////////////////////////////////////////////////////////////////
function is_internal_point(lar::Lar, ray_origin::PointNd, ray_dir::PointNd)
  pdim=size(lar.V,1)
  num_intersections=0
  for I in 1:length(lar.C[ pdim==2 ? :EV : :FV])
    hit,t=ray_face_intersection(ray_origin, ray_dir, lar, I)
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
  # VIEWCOMPLEX(lar)

  pdim=size(lar.V, 1)
  @assert(pdim==2 || pdim==3)

	# TODO: reintroduce space index to avoid O(N^2) complexity

  b1,b2=lar_bounding_box(lar; only_used_vertices=true)
  
  # i want to move out of the complex
  move_out = 3.0 * LinearAlgebra.norm(b2 - b1)   

  @show(b1,b2)
  @show(move_out)

  for attempt in 1:max_attempts

    # choose a random point inside the box and a random direction
    inside_bbox = [random_float(b1[I],b2[I]) for I in 1:pdim ]

    external_point = inside_bbox + move_out * random_dir(pdim)
    ray_dir    = normalized(inside_bbox-external_point)

    @show("Trying", external_point, ray_dir)

    distances=PointNd()

    for I in 1:length(lar.C[ pdim==2 ? :EV : :FV])
      hit,distance=ray_face_intersection(external_point, ray_dir, lar, I)
      if !isnothing(hit)
        push!(distances,distance)
        @show("hit", external_point,ray_dir,distance)
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
      @show(distances)
      best_delta,internal_point=nothing,nothing
      for I in 1:length(distances)-1
        delta=distances[I+1]-distances[I]
        if ((I % 2) ==1) && delta>LAR_DEFAULT_ERR # only if is inside and delta is not too small
          distance=distances[I]+delta/2
          if isnothing(best_delta) || delta>best_delta
            internal_point=external_point+ray_dir*distance
            @show("new best", external_point, ray_dir, distance,internal_point, best_delta,delta)
            best_delta=delta
          else
            @show("not the best",delta)
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
    @show(internal_point, ray_dir, distances)
    return internal_point, ray_dir, distances

  end

	error("cannot find internal point")
end

# /////////////////////////////////////////////////////////////////////
function Union(v::Vector{Bool})::Bool       
  return any(v) 
end

function Intersection(v::Vector{Bool})::Bool
  return all(v) 
end

function Difference(v::Vector{Bool})::Bool 
  return v[1] && !any(v[2:end]) 
end

function Xor(v::Vector{Bool})::Bool  
  return (length([it for it in v if it]) % 2)==1  
end

export Union, Intersection, Difference, Xor

# /////////////////////////////////////////////////////////////////////
function SELECT_ATOMS(lar::Lar, sel::Cell)::Lar

  pdim=size(lar.V, 1)

  if pdim==2

    Fsel, Esel=Set{Int}(), Set{Int}()
    
    for F in sel
      push!(Fsel,F)
      for E in lar.C[:FE][F]
        push!(Esel,E)
      end
    end

    ret=Lar(lar.V, Dict(
        :FE => Cells(),
        :EV => Cells())
    )

    ret.mapping=Dict(
      :F => Dict{Int,Int}(), 
      :E => Dict{Int,Int}()
    )

    # add edges
    Emap=Dict{Int,Int}()
    for Eold in Esel
      push!(ret.C[:EV], [ v_index for v_index in lar.C[:EV][Eold] ])
      Enew=length(ret.C[:EV])
      ret.mapping[:E][Enew]=Eold
      Emap[Eold]=Enew
    end

    # add faces
    Fmap=Dict{Int,Int}()
    for Fold in Fsel
      push!(ret.C[:FE], [ Emap[Eold] for Eold in lar.C[:FE][Fold] ])
      Fnew=length(ret.C[:FE])
      ret.mapping[:F][Fnew]=Fold
      Fmap[Fold]=Fnew
    end

    ret.C[:FV]=compute_FV(ret)
    return ret

  else

    Csel, Fsel, Esel=Set{Int}(),Set{Int}(), Set{Int}()
    for C in sel
      push!(Csel,C)
      for F in lar.C[:CF][C]
        push!(Fsel,F)
        for E in lar.C[:FE][F]
          push!(Esel,E)
        end
      end
    end

    ret=Lar(
      lar.V, 
      Dict(
        :CF => Cells(), 
        :FE => Cells(),
        :EV => Cells())
    )

    ret.mapping=Dict(
      :C => Dict{Int,Int}(),
      :F => Dict{Int,Int}(), 
      :E => Dict{Int,Int}()
    )

    # add edges
    Emap=Dict{Int,Int}()
    for Eold in Esel
      push!(ret.C[:EV], [ v_index for v_index in lar.C[:EV][Eold] ])
      Enew=length(ret.C[:EV])
      ret.mapping[:E][Enew]=Eold
      Emap[Eold]=Enew
    end

    # add faces
    Fmap=Dict{Int,Int}()
    for Fold in Fsel
      push!(ret.C[:FE], [ Emap[Eold] for Eold in lar.C[:FE][Fold] ])
      Fnew=length(ret.C[:FE])
      ret.mapping[:F][Fnew]=Fold
      Fmap[Fold]=Fnew
    end

    # add atoms
    for Cold in Csel
      push!(ret.C[:CF], [ Fmap[Fold] for Fold in lar.C[:CF][Cold] ])
      Cnew=length(ret.C[:CF])
      ret.mapping[:C][Cnew]=Cold
    end

    ret.C[:FV]=compute_FV(ret)
    ret.C[:CV]=compute_CV(ret)

    return ret
  end

end

# //////////////////////////////////////////////////////////////////////////////
function ATOMS(lar::Lar)::Vector{Lar}
  pdim=size(lar.V,1)
  return Vector{Lar}([SELECT_ATOMS(lar, [A]) for A in eachindex(lar.C[pdim==2 ? :FE : :CF])])
end
export ATOMS

# ////////////////////////////////////////////////////////////////
function BOOL(arrangement::Lar; input_args=[], bool_op=Union, debug_mode=false)::Lar
  
  # if you want to see atoms set debug_mode=true
  atoms=ATOMS(arrangement)

  internal_points=[find_internal_point(atom) for atom in atoms] 
  
  bool_matrix=zeros(Bool,length(atoms),length(input_args))
  
  for (A,(atom, (internal_point, ray_dir, distances))) in enumerate(zip(atoms, internal_points))
    # @show(internal_point, ray_dir, distances)
  
    # find for each input arg if the internal point is inside or outside
    for (I,input_arg) in enumerate(input_args)
      bool_matrix[A,I]=is_internal_point(input_arg, internal_point, ray_dir)
    end
  
    # show internal / external points for classification
    if debug_mode

      # todo 2d
      @assert(pdim==3)

      viewer=Viewer()

      # points
      begin
        points, colors = Vector{Float32}(),Vector{Float32}()
        append!(colors,RED)
        append!(points, internal_point)
        for distance in distances
          append!(colors, GREEN )
          append!(points, internal_point+distance*ray_dir)
        end
        render_points(viewer, points, colors=colors, point_size=8)
      end

      # lines
      begin
        lines , colors = Vector{Float32}(),Vector{Float32}()
        for distance in distances
          append!(colors,  YELLOW); append!(lines,  internal_point)
          append!(colors,  YELLOW); append!(lines,  internal_point+distance*ray_dir)
        end
        render_lines(viewer, lines,  colors=colors, line_width=DEFAULT_LINE_WIDTH)
      end

      explanation=join([string(it) for it in bool_matrix[A,:]], " ")
      # render_text(viewer, explanation, center=internal_point, color=BLACK)

      VIEWCOMPLEX(viewer, atom, show=["V", "EV"], explode=[1.0,1.0,1.0], title=explanation)
    end
  end
  
  sel=Cell()
  for (A,row) in enumerate(eachrow(bool_matrix))
    if bool_op(collect(row))
      append!(sel, A)
    end
  end

  return SELECT_ATOMS(arrangement, sel)

end

export BOOL








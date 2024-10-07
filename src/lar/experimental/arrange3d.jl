

# //////////////////////////////////////////////////////////////
function arrange3d_v2(lar::Lar; debug_mode=true)

  POINTS3D, EV, FE=PointsDB(), Cells(),Cells()

  num_faces=length(lar.C[:FV])
  for (F1, fv) in enumerate(lar.C[:FV])
    println("# // fragmenting face=$(F1) / $(num_faces)")

    points2d, segments = PointsDB(), Cells()

    projector, roi_3d, roi_2d=nothing, nothing, nothing
    F2s=collect(1:length(lar.C[:FV]))
    delete!(F2s,F1)
    for (K,F) in enumerate([F1 ; F2s ])
  
      fv=lar.C[:FV][F]
      face_points3d=lar.V[:,fv]
      face_box3d=bbox_create(face_points3d)
  
      if K==1
        plane=plane_create(face_points3d)
        projector=project_points3d(face_points3d, double_check=true) # scrgiorgio: remove double check 
        roi_3d=face_box3d
      else
        if !bbox_intersect(roi_3d, face_box3d) 
          continue 
        end
      end
  
      face_points_2d=Vector{Vector{Float64}}()
      for E in lar.C[:FE][F]
        a,b=lar.C[:EV][E]
        p1,p2=[Vector{Float64}(it) for it in eachcol(projector(lar.V[:,[a,b]]))]
        append!(face_points_2d,[p1, p2])
      end

      if K==1
        roi_2d=bbox_create(face_points_2d)
      end
    
      for (p1,p2) in face_points_2d
        if K==1 || !bbox_intersect(roi_2d, bbox_create([p1,p2]))
          a=add_point(points2d, p1)
          b=add_point(points2d, p2)
          push!(segments,[a,b]) 
        end
      end
  
    end
  
    segments=[it for it in simplify_cells(segments) if length(it)==2]
    tin = Triangulate.TriangulateIO()
    tin.pointlist  = ToPoints(points2d)
    tin.segmentlist = hcat(segments...)  # constrained triangulation
  
    if debug_mode
      VIEWEDGES(tin.pointlist, tin.segmentlist, title="arrange3d_v2 / Triangulate input")
    end
  
    # https://github.com/JuliaGeometry/Triangulate.jl/blob/8e0fc2c0fb58ffb3ae351afcca5da375522d2e84/docs/src/triangle-h.md?plain=1#L193
    # -p Triangulates a Planar Straight Line Graph, 
    # -Q for quiet
    tout=nothing
    while true
      try
        tout, __ = Triangulate.triangulate("pQ", tin) 
        break
      catch TriangulateError
        println("WARNING Triangulate.triangulate failed, so perturnbing the points")
        tin.pointlist=ToPoints([p + rand(2) * LAR_ARRANGE2D_PERTURBATION for p in eachcol(lar.V)])
      end
    end
  
    if debug_mode
      VIEWTRIANGLES(tout.pointlist, tout.trianglelist, title="arrange3d_v2 / Triangulate output")
    end
  
    segments  = Set(simplify_cells([Cell([a,b])   for (a,b  ) in eachcol(tout.segmentlist)]))
    triangles =     simplify_cells([Cell([a,b,c]) for (a,b,c) in eachcol(tout.trianglelist)])
  
    # unproject and round vertices (in 3D! so that faces can glue together)
    begin
      vmap=Dict()
      unprojected=projector(tout.pointlist, inverse=true)
      for P in 1:size(tout.pointlist,2)
        p3d=unprojected[:,P]
        r3d=round_vector(Vector{Float64}(p3d), digits=LAR_ARRANGE2D_ROUND)
        vmap[P]=add_point(POINTS3D, r3d)
      end
    end
  
    # build faces (vertex indices are now referring to rounded 3d)
    begin
      face_EV=Cells()
      face_FE=Cells()
      for cycle in find_triangles_cycles(triangles, segments)
        fe=Cell()
        for (C,(a,b)) in enumerate(cycle)
          a,b= (C==length(cycle)) ? cycle[end] : [cycle[C][1],cycle[C+1][1]]
          a,b=vmap[a],vmap[b]
          push!(face_EV, [a,b]) # will simplify later
          push!(fe, length(face_EV))
        end
        push!(face_FE, fe)
      end
      face_EV=[it for it in simplify_cells(face_EV) if length(it)>=2]
      face_FE=[it for it in simplify_cells(face_FE) if length(it)>=3]
      append!(EV, face_EV)
      append!(FE, face_FE)
    end
  end

  POINTS3D=get_points(POINTS3D)
  lar=Lar(POINTS3D, Dict( :EV => EV, :FE => FE))

  # if you had a problem with an edge, you should be able to debug here
  # debug_edge(lar,21,26)

  lar=SIMPLIFY(lar)
  COMPUTE(lar,:FV)

  if debug_mode
    VIEWCOMPLEX(lar, explode=[1.0,1.0,1.0], show=["V","FV","Vtext"], title="arrange3d_v2 / all faces fragmented")
  end  

  # find atoms
  lar.C[:CF]=lar_find_atoms(lar.V, cycles, debug_mode=debug_mode)

  # show atoms
  if debug_mode
    VIEWCOMPLEX(lar, show=["FV","atom"], explode=[1.2,1.2,1.2], title="arrange3d_v2 / atoms")
  end

  return lar
end
export arrange3d_v2



# ///////////////////////////////////////////////////////////////////////////
# Newell's method, works for concave too and it is oriented
function compute_oriented_newell_normal(loop::Vector{Vector{Float64}})::Vector{Float64}
    n = [0.0, 0.0, 0.0]
    for (I, (x1, y1, z1)) in enumerate(loop)
        x2, y2, z2 = loop[mod1(I+1, length(loop))]
        n[1] += (y1 - y2) * (z1 + z2)
        n[2] += (z1 - z2) * (x1 + x2)
        n[3] += (x1 - x2) * (y1 + y2)
    end
    return n / LinearAlgebra.norm(n)
end

# ///////////////////////////////////////////////////////////////////////////
# tgw angle computation
function compute_oriented_angle(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64})::Float64
    a = a / LinearAlgebra.norm(a)
    b = b / LinearAlgebra.norm(b)
    angle = atan(dot(cross(a, b), c), dot(a, b))
    return rad2deg(angle)
end

# ///////////////////////////////////////////////////////////////////////////
function test_compute_normal()
    @test compute_oriented_newell_normal([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]) ≈ [0.0, 0.0, +1.0]
    @test compute_oriented_newell_normal([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]]) ≈ [0.0, 0.0, -1.0]
end

function test_get_oriented_angle()
    a, edge = [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]
    @test compute_oriented_angle(a, a, edge) ≈ 0.0
    @test compute_oriented_angle(a, [.00, -1.0,  0.0], edge) ≈ -90.0
    @test compute_oriented_angle(a, [.00, -1.0, -1.0], edge) < -90.0
    @test compute_oriented_angle(a, [.00, +1.0,  0.0], edge) ≈ +90.0
    @test compute_oriented_angle(a, [.00, +1.0, -1.0], edge) > +90.0

    a, edge = [0.0, 0.0, +1.0], [+1.0, 0.0, 0.0]
    @test compute_oriented_angle(a, a, edge) ≈ 0.0
    @test compute_oriented_angle(a, [0.0, +1.0,  0.0], edge) ≈ -90.0
    @test compute_oriented_angle(a, [0.0, +1.0, -1.0], edge) < -90.0
    @test compute_oriented_angle(a, [0.0, -1.0,  0.0], edge) ≈ +90.0
    @test compute_oriented_angle(a, [0.0, -1.0, -1.0], edge) > +90.0

    a, edge = [0.0, 0.0, -1.0], [+1.0, 0.0, 0.0]
    @test compute_oriented_angle(a, a, edge) ≈ 0.0
    @test compute_oriented_angle(a, [0.0, -1.0,  0.0], edge) ≈ -90.0
    @test compute_oriented_angle(a, [0.0, -1.0, +1.0], edge) < -90.0
    @test compute_oriented_angle(a, [0.0, +1.0,  0.0], edge) ≈ +90.0
    @test compute_oriented_angle(a, [0.0, +1.0, +1.0], edge) > +90.0

    a, edge = [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]
    @test compute_oriented_angle(a, a, edge) ≈ 0.0
    @test compute_oriented_angle(a, [0.0,  1.0,  0.0], edge) ≈ -90.0
    @test compute_oriented_angle(a, [0.0,  1.0,  1.0], edge) < -90.0
    @test compute_oriented_angle(a, [0.0, -1.0,  0.0], edge) ≈ +90.0
    @test compute_oriented_angle(a, [0.0, -1.0, +1.0], edge) > +90.0
end


# ////////////////////////////////////////////////////////////////////////
function lar_connected_components(seeds::Vector{Int}, get_connected::Function)::Vector{Vector{Int}}
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

# ////////////////////////////////////////////////////////////////////////
function lar_find_atoms(V::Points, cycles::Cycles; debug_mode=false)::Cells

  if true
    println("Cycles")
    for (C,cycle) in enumerate(cycles)
      println(cycle, " # ",C)
    end
  end

  #0.5 1.0 0.0;# 187
  #0.0 0.5 -0.5;# 224

  num_cycles=length(cycles)

  connections::Dict{Int,Dict} = Dict(
      ( A => Dict() for A in 1:num_cycles)...,
      (-A => Dict() for A in 1:num_cycles)...
  )

  @assert(all([abs(k)>=1 && abs(k)<=num_cycles for k in keys(connections)]))
  for (A, Acycle) in enumerate(cycles)
    @assert(A>=1 && A<=num_cycles)
    for (a, b) in Acycle
      adjacent1, adjacent2 = Vector{Any}(), Vector{Any}()
      
      # other faces incidend to the same edge
      Bs=collect(unique([B for (B, Bcycle) in enumerate(cycles) if (B!=A) && (([a,b] in Bcycle) || ([b,a] in Bcycle))]))
      if length(Bs)==0
        println("PROBLEM WITH edge ",a," ",b)
        @assert(false)
      end

      for B in Bs

        @assert(B>=1 && B<=num_cycles)

        # coherent loop
        if [b, a] in cycles[B]
            Bcycle, B = cycles[B], +B
        elseif [a, b] in cycles[B]
            Bcycle, B = reverse_cycle(cycles[B]), -B
        else
            @assert false
        end
        An = compute_oriented_newell_normal([V[:,first] for (first, ___) in Acycle])
        Bn = compute_oriented_newell_normal([V[:,first] for (first, ___) in Bcycle])
        Ev = V[:,b] - V[:,a]
        angle1=compute_oriented_angle(+1.0 .* An, +1.0 .* Bn, +1.0 .* Ev)
        angle2=compute_oriented_angle(-1.0 .* An, -1.0 .* Bn, -1.0 .* Ev)
        push!(adjacent1, (face=+B, angle=angle1))
        push!(adjacent2, (face=-B, angle=angle2))
      end

      adjacent1=collect(sort(adjacent1,by=value->value.angle));@assert(length(adjacent1)>0)
      adjacent2=collect(sort(adjacent2,by=value->value.angle));@assert(length(adjacent2)>0)
      ev=collect(sort([a,b]))
      connections[+A][ev]=adjacent1  
      connections[-A][ev]=adjacent2

    end
  end

  # print connections for debugging
  if debug_mode
    for (F,  adjacent_per_edge) in connections
      @assert(abs(F)>=1 && abs(F)<=num_cycles)
      # println("F=", F)
      for (ev, adjacents) in adjacent_per_edge
        # println("  ev=",ev, adjacents)
      end
    end
  end

  # topology checks on connections
  begin
    for (A,  adjacent_per_edge) in connections
      for (ev, adjacents) in adjacent_per_edge
        for adj in adjacents
          B=adj.face
          @assert B in [it.face for it in connections[A][ev]]
          @assert A in [it.face for it in connections[B][ev]]
        end
      end
    end
  end

  # find best connections by angles (TGW)
  begin
    best_connections=Dict{Int, Vector{Int}}()
    for (A,  adjacent_per_edge) in connections
      best_connections[A]=remove_duplicates([adjacents[1].face for (ev, adjacents) in adjacent_per_edge])
    end
  end

  @assert(all([abs(F)>=1 && abs(F)<=num_cycles for F in keys(best_connections)]))
  atoms=lar_connected_components(collect(keys(best_connections)), cur -> best_connections[cur])

  # atoms topology checks
  begin

    @assert(all([abs(F)>=1 && abs(F)<=num_cycles for F in vcat(atoms...)]))

    for F in keys(connections)

      # all index faces should be fine

      # one signed face should be in only one atom
      @assert(length([atom for atom in atoms if F in atom])==1)

      # if there is +F in one atom, there cannot be -F
      @assert(length([atom for atom in atoms if F in atom && -F in atom])==0)

    end
  end

  # I do not need the sign anymore (could be useful for debugging)
  atoms=[collect(sort([abs(jt) for jt in it])) for it in atoms]


  # @show(atoms)

  # topology check
  begin
    num_full_cell_per_face= Dict{Int,Int}()
    for (A, atom) in enumerate(atoms)

      # each atom must not contain the same face twice
      if length(Set(atom))!=length(atom)
        @assert(false)
      end

      for F in atom
        if !haskey(num_full_cell_per_face,F) num_full_cell_per_face[F]=0 end
        num_full_cell_per_face[F]+=1
      end
    end

    # @show(num_full_cell_per_face)

    # no hanging face, and a face cannot be shared by more than 2-full-dim cells
    # note: consider the outer atom too
    @assert(all([v==2 for (k,v) in num_full_cell_per_face]))

  end

  return atoms
end

# //////////////////////////////////////////////////////////////
function debug_edge(lar::Lar, a,b; enlarge=nothing)

  @show(lar)
  for (I,fv) in enumerate(lar.C[:FV])
    println("fv ", fv, " # ",I)
  end
  
  selected_vertices=[a,b]

  # take any vertex near by
  if !isnothing(enlarge)
    for (P,p) in enumerate(eachcol(lar.V))
      if LinearAlgebra.norm(p - lar.V[:,a])<enlarge || LinearAlgebra.norm(p - lar.V[:,b])<enlarge push!(selected_vertices,P) end
    end
  end

  # show all faces touching the vertices
  selected_faces=Cell()
  for (F,fv) in enumerate(lar.C[:FV])
    inters=intersect(Set(selected_vertices),Set(fv))
    if length(inters)>0
      push!(selected_faces,F)
    end
  end

  selected_faces=normalize_cell(selected_faces)
  VIEWCOMPLEX(SELECT(lar,selected_faces), show=["FV", "Vtext"], explode=[1.0,1.0,1.0], face_color=TRANSPARENT, title="debug edge")

end


# //////////////////////////////////////////////////////////////////////////////
function guess_boundary_faces(lar::Lar, faces::Vector; max_attempts=1000)::Vector{Int}
  for attempt in 1:max_attempts
    b1,b2=lar_bounding_box(lar; only_used_vertices=true)
    move_out = 3*LinearAlgebra.norm(b2 - b1)   
    inside_bbox = [random_float(b1[I],b2[I]) for I in 1:3 ]
    external_point = inside_bbox + move_out * random_dir()
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
function compute_atom_bbox(lar::Lar, atom::Cell)
  points=Vector{Vector{Float64}}()
  for F in atom
    fv=lar.C[:FV][F]
    append!(points,[p for p in eachcol(lar.V[:,fv])])
  end
  return bbox_create(points)
end

# ////////////////////////////////////////////////////////////////
function arrange3d_v2_split(lar::Lar)::Tuple{Lar,Lar}

  atoms=[cf for cf in lar.C[:CF]]

  # connected atoms
  begin
    # @show(atoms)
    components=lar_connected_components(collect(eachindex(atoms)), A -> [B for B in eachindex(atoms) if A!=B && length(intersect(Set(atoms[A]),Set(atoms[B])))>0 ])
    components=[ remove_duplicates([ atoms[idx] for idx in component]) for component in components]
    # NOTE: components are normalized (no repetition and sorted)
  end

  # topology check on components
  begin
    for component in components
      @assert(length(remove_duplicates(component))==length(component))
    end
  end

  lar_outers=lar_copy(lar);lar_outers.C[:CF]=[]
  lar_inners=lar_copy(lar);lar_inners.C[:CF]=[]
  
  for (C,component) in enumerate(components)

    faces=remove_duplicates(vcat(component...))
    num_atoms=length(component)
    # println("# components ",C, " num_atoms=", num_atoms, " faces=",faces)
    # println(component)
    @assert(length(component)>=2)

    outer,inners=nothing,[]

    # there is one outer cell, and one inner cell and they must be the same
    if num_atoms==2
      @assert(component[1]==component[2])
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
function arrange3d_v2_inners(lar::Lar)::Lar
  outers,inners = arrange3d_v2_split(lar)
  return inners
end

# ////////////////////////////////////////////////////////////////
function arrange3d_v2_outers(lar::Lar)::Lar
  outers,inners=arrange3d_v2_split(lar)
  return outers
end







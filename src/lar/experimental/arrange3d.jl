# //////////////////////////////////////////////////////////////
function arrange3d_v2(src::Lar; debug_mode=true, debug_face=:none, debug_edge=nothing)

  dst=Lar()
  dst.C[:EV]=Cells()
  dst.C[:FE]=Cells()
  points3d=PointsDB()

  num_faces=length(src.C[:FV])
  for (F1, fv) in enumerate(src.C[:FV])
    println("# fragmenting face=$(F1) / $(num_faces)")

    points2d = PointsDB()

    # F1
    begin
      F1_points_3d=src.V[:,src.C[:FV][F1]]
      plane=plane_create(F1_points_3d)
      projector=project_points3d(F1_points_3d, double_check=true) # scrgiorgio: remove double check 
      F1_box_3d = bbox_create(F1_points_3d)
      F1_segments = Cells()
      for E in src.C[:FE][F1]
        a,b=src.C[:EV][E]
        p1,p2=[PointNd(it) for it in eachcol(projector(src.V[:,[a,b]]))]
        a=add_point(points2d, p1)
        b=add_point(points2d, p2)
        push!(F1_segments, [a, b])
      end
      F1_box_2d        = bbox_create(get_points(points2d))
      F1_points_2d_row = BYROW(get_points(points2d))
      F1_is_inside = p -> (classify_point(p, F1_points_2d_row, F1_segments)  != "p_out")
    end

    # F2s
    begin
      F2_segments=Cells()
      for F2 in collect(1:num_faces)

        if F2==F1
          continue
        end

         # quick discard not intersecting faces
        F2_box_3d=bbox_create(src.V[:,src.C[:FV][F2]])
        if !bbox_intersect(F1_box_3d,F2_box_3d) 
          continue 
        end


        # find intersection on the main face
        begin
          hits=PointsNd()
          for E2 in src.C[:FE][F2]
            a,b = src.C[:EV][E2]
            p1,p2=[PointNd(it) for it in eachcol(src.V[:,[a,b]])]
            hit, t = plane_ray_intersection(p1, normalized(p2-p1), plane)
            if isnothing(hit) continue end
            @assert(t>0) # must cross the plane (be on opposite sides of the plane)
            if t<LinearAlgebra.norm(p2-p1)
              push!(hits, projector(hcat(hit))[:,1])
            end
          end
        end

        # scrgiorgio: WRONG. I should take all hits, start from the longest distance and build segment on,off,on,off
        if length(hits)>0

          # not sure what this means... should always F1 cross F1 or can just touch?
          if length(hits)==1
            continue
          end

          # TODO: ordering here

          for (a,b) in zip(hits[1:end-1],hits[2:end])
            if bbox_intersect(F1_box_2d, bbox_create([a,b])) 
              A=add_point(points2d, a)
              B=add_point(points2d, b)
              push!(F2_segments, [A,B]) 
            end
          end
        end

      end
    end
  
    # triangulate
    begin
      segments=[it for it in simplify_cells([F1_segments ;  F2_segments]) if length(it)==2] # remove degenerate (probably not needed since I am not rounding)
      tin = Triangulate.TriangulateIO()
      tin.pointlist  = get_points(points2d)
      tin.segmentlist = hcat(segments...)  # constrained triangulation
    
      if debug_face==F1 || debug_face==:all
        VIEWEDGES(tin.pointlist, tin.segmentlist, title="arrange3d / 2d triangulate input face=$(F1)")
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
          println("WARNING Triangulate.triangulate failed, so perturbing the points")
          tin.pointlist=ToPoints([p + rand(2) * LAR_EXPERIMENTAL_ARRANGE_PERTURBATION for p in eachcol(get_points(points2d))])
        end
      end
    
      if debug_face==F1 || debug_face==:all
        VIEWEDGES(tout.pointlist, tout.segmentlist, title="arrange3d / 2d triangulate output face=$(F1)", explode=[1.2,1.2,1.2])
        #VIEWEDGES(tout.pointlist, tout.trianglelist, title="arrange3d / 2d triangulate output face=$(F1)")
      end
    end

    # round vertices so that faces can glue together)
    begin
      v_index_2d_to_3d=Dict()
      v_index_3d_to_2d=Dict()
      unprojected=projector(tout.pointlist, inverse=true)
      for index_2d in 1:size(tout.pointlist,2)
        p3d=unprojected[:,index_2d]
        r3d=round_vector(PointNd(p3d), digits=LAR_EXPERIMENTAL_ARRANGE_ROUND)
        index_3d=add_point(points3d, r3d)
        v_index_2d_to_3d[index_2d]=index_3d
        v_index_3d_to_2d[index_3d]=index_2d
      end
      dst.V=get_points(points3d)
    end

    # find cycles (in the original triangulate space, i.e. 2d)
    begin
      cycles=Cycles()

      V             = tout.pointlist
      segments      = simplify_cells([Cell([a,b])   for (a,b  ) in eachcol(tout.segmentlist)])
      triangles     = simplify_cells([Cell([a,b,c]) for (a,b,c) in eachcol(tout.trianglelist)])
      #VIEWTRIANGLES(V, triangles, title="triangle F=$(F1)", explode=[2.2,2.2,2.2])


      # do the rounding, need to do here because no one can guarantee I have cycles above
      # (instead we are hoping cycles will be correct in the rounded space)
      begin

        triangles = [[v_index_3d_to_2d[v_index_2d_to_3d[index]] for index in it]  for it in triangles ]
        segments  = [[v_index_3d_to_2d[v_index_2d_to_3d[index]] for index in it]  for it in segments  ]
          
        # remove degenerate due to the rounding
        triangles = [it for it in simplify_cells(triangles) if length(it)==3] 
        segments = Set([it for it in simplify_cells(segments) if length(it)==2])
      end

      # if something anything could go wrong here, and not sure what to do
      # i am betting that the rounding keeps some topological persistency
      begin
        adjacent_triangles=find_adjacents_cells(triangles, 2, segments)
        groups=find_groups_of_cells(adjacent_triangles)
      end

      for (A, triangle_ids) in enumerate(groups)

        # each group will form a face (can be holed and non-convex)
        inside_area,outside_area=0.0,0.0

        complex_face=Cells()
        for triangle_id in triangle_ids 
          u,v,w = triangles[triangle_id]

          # need to exclude outside F1
          centroid=(V[:,u] + V[:,v] + V[:,w])/3.0
          area=GetTriangleArea(V[:,u],V[:,v], V[:,w])
          if F1_is_inside(centroid) 
            inside_area+=area
          else 
            outside_area+=area 
          end

          for (a,b) in [ [u,v], [v,w], [w,u] ]
            a,b = normalize_cell([a,b])
            if [a,b] in segments
              push!(complex_face, [a,b])
            end
          end
        end

        # if you want to debug the group
        #if debug_face==F1
        #  VIEWTRIANGLES(V, Cells([triangles[tt] for tt in triangle_ids]), title="triangle group $(F1)")
        #end

        # println("   inside_area=$(inside_area) outside_area=$(outside_area)")
        if inside_area>0 

          if outside_area>0
            println("WARNING ambiguity F=$(F1) inside_area=$(inside_area) outside_area=$(outside_area)")
            # VIEWTRIANGLES(V, Cells([triangles[tt] for tt in triangle_ids]), title="xxx $(F1)")
          end

          if inside_area> outside_area
            complex_face=simplify_cells(complex_face)
            for cycle in find_vcycles(complex_face)

              # go to the rounded world, so something can get filtered
              cycle=normalize_cycle( [[v_index_2d_to_3d[a],v_index_2d_to_3d[b]] for (a,b) in cycle]) 
              if length(cycle)>=3
                push!(cycles, cycle)
              end
            end
          end

        end

      end
    end
  
    # build faces (vertex indices are now referring to rounded 3d)
    begin
      sel=Cell()
      for cycle in cycles
        fe=Cell()
        for (C,(a,b)) in enumerate(cycle)
          a,b= (C==length(cycle)) ? cycle[end] : [cycle[C][1],cycle[C+1][1]]
          push!(dst.C[:EV], [a,b]) # will simplify later
          push!(fe, length(dst.C[:EV]))
        end
        push!(dst.C[:FE], fe)
        push!(sel, length(dst.C[:FE]))
      end
    end

    vids=vcat([[a for (a,b) in cycle] for cycle in cycles]...)
    if debug_face==F1 || debug_face==:all || (!isnothing(debug_edge) && debug_edge[1] in vcat(vids...) && debug_edge[2] in vids)
      VIEWCOMPLEX(SELECT(dst, sel), show=["V","FV","Vtext"], title="arrange3d / 3d face face=$(F1)")
    end  

  end

  # remove spurios/isolated vertex on cell boundaries (due to the triangulate library)
  begin

    VE=compute_VE(dst)
    for (v_index,ve) in enumerate(VE)
      # only 3 vertices, it should be >=3
      if length(ve)==2
        println("Removing isolated vertex $(v_index)")
        E1,E2=ve
        (a,b),(c,d)= dst.C[:EV][E1], dst.C[:EV][E2]
        a,b = [it for it in [a,b,c,d] if it!=v_index]
        @assert(a!=b && a!=v_index && b!=v_index)
        dst.C[:EV][E1]=[a,b]
        dst.C[:EV][E2]=[a,b] # replicate so it will be removed by SIMPLIFY below
      end
    end

  end

  dst=SIMPLIFY(dst)

  dst.C[:FV]=compute_FV(dst)

  # if you had a problem with an edge, you should be able to debug here
  if !isnothing(debug_edge)
    VIEWEDGE(dst, debug_edge)
  end

  if debug_mode
    VIEWCOMPLEX(dst, explode=[1.0,1.0,1.0], show=["V","EV","Vtext"], title="arrange3d / 3d ALL faces")
  end  

  # find atoms
  begin
    cycles=Cycles()
    for fe in dst.C[:FE]
      cell_EV=Cells([dst.C[:EV][E] for E in fe])
      append!(cycles, find_vcycles(cell_EV))
    end
    dst.C[:CF]=lar_find_atoms(dst.V, cycles, debug_mode=debug_mode)
  end

  # show atoms
  if debug_mode
    VIEWCOMPLEX(dst, show=["CV"], explode=[1.2,1.2,1.2], title="arrange3d / ALL atoms")
  end

  return dst
end
export arrange3d_v2



# ///////////////////////////////////////////////////////////////////////////
# Newell's method, works for concave too and it is oriented
function compute_oriented_newell_normal(loop::AbstractPointsNd)::PointNd
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
function compute_oriented_angle(a::PointNd, b::PointNd, c::PointNd)::Float64
    a = a / LinearAlgebra.norm(a)
    b = b / LinearAlgebra.norm(b)
    angle = atan(dot(cross(a, b), c), dot(a, b))

    ret=rad2deg(angle)

    # I think it will not be a problem if the first face I see is 180 or -180 it means they are two flat faces connected each others
    # @assert(ret!=180.0 || ret==-180.0) # ambiguity otherwise
    return ret
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

# ////////////////////////////////////////////////////////////////////////
function lar_find_atoms(V::Points, cycles::Cycles; debug_mode=false)::Cells

  #0.5 1.0 0.0;# 187
  #0.0 0.5 -0.5;# 224

  num_cycles=length(cycles)

  connections::Dict{Int,Dict} = Dict(
      ( A => Dict() for A in 1:num_cycles)...,
      (-A => Dict() for A in 1:num_cycles)...
  )

  for (A, Acycle) in enumerate(cycles)
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
      #println("F=", F)
      for (ev, adjacents) in adjacent_per_edge
        #println("  ev=",ev, adjacents)
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
    best_connections=Dict{Int, Cell}()
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


# //////////////////////////////////////////////////////////////////////////////
function guess_boundary_faces(lar::Lar, faces::Vector; max_attempts=1000)::Cell
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
  points=PointsNd()
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
function arrange3d_v2_inners(lar::Lar)::Lar
  outers,inners = arrange3d_v2_split(lar)
  return inners
end

# ////////////////////////////////////////////////////////////////
function arrange3d_v2_outers(lar::Lar)::Lar
  outers,inners=arrange3d_v2_split(lar)
  return outers
end







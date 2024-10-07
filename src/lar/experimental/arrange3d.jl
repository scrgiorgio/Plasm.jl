

# /////////////////////////////////////////////////////////////////////
function transform_triangles(V::Points, in_segments::Cells, in_triangles::Cells; 
  segment_map::Dict{Cell,Cells}=Dict{Cell,Cells}(), triangle_map=Dict(),vertex_map=Dict())

  out_segments=Cells()
  segment_map=Dict(normalize_cell(k) => v  for (k,v) in segment_map)
  for Old::Cell in in_segments
    for New::Cell in get(segment_map, normalize_cell(Old), [Old])
      New=normalize_cell(New)
      # apply vertex map
      New=[get(vertex_map,v_index,v_index) for v_index in New] 
      push!(out_segments, New)
    end
  end
  out_segments::Cells=[it for it in simplify_cells(out_segments) if length(it)==2]

  out_triangles=Cells()
  triangle_map=Dict(normalize_cell(k) => v  for (k,v) in triangle_map )
  for Old::Cell in in_triangles
    for New::Cell in get(triangle_map, normalize_cell(Old), [Old])
      New=normalize_cell(New)
      # apply vertex map
      New=[get(vertex_map,v_index,v_index)  for v_index in New]
      push!(out_triangles, New)
    end
  end 
  out_triangles::Cells=[it for it in simplify_cells(out_triangles) if length(it)==3]

  return V, out_segments, out_triangles
end

# /////////////////////////////////////////////////////////////////////
function remove_small_or_skewed_triangles(V, segments::Cells, triangles::Cells; epsilon=LAR_ARRANGE2D_SMALL_TRIANGLES_ERR)

  # go in order by number of edges <= epsilon
  for required_num_small in [3,2,1,0]

    for (a,b,c) in triangles

      triangle_info=GetTriangleInfo(V, a,b,c)
      edge_lengths=triangle_info[:edge_lengths]
      num_small=length([it for it in edge_lengths if it<epsilon])

      # I want to go in order from the smalled stuff to the largest stuff
      if num_small!=required_num_small
        continue
      end

      # [OK]
      # collapsing triangle and its edgesin a point
      # **hoping (ab),(bc),(c,d) are not `boundaries` because they disappear**
      # NOTE: this will remove the 3 edges and 3 triangles since they will become degenerate
      if num_small==3
        println("Removing triangle because too small overall triangle=", [a,b,c], " edge_lengths=", edge_lengths)
        # VIEWTRIANGLES(V, triangles, explode=[1.8,1.8,1.8])
        return transform_triangles(V, segments, triangles, vertex_map=Dict( b=> a, c=>a )) 
      end

      # [OK]
      # collapsing triangle and its edges in a point
      # 3rd edge has module in the range [epsilon,2epsilon] and so I think it can be removed as well
      # **hoping (ab),(bc),(c,d) are not `boundaries` because they disappear**
      # NOTE: this will remove the 3 edges and 3 triangles since they will become degenerate
      if num_small==2
        println("Removing triangle because 2 small edges and 1 slightly more than small=", [a,b,c], " edge_lengths=", edge_lengths)
        # VIEWTRIANGLES(V, triangles, explode=[1.8,1.8,1.8])
        return transform_triangles(V, segments, triangles, vertex_map=Dict( b=> a, c=>a )) 
      end

      # [OK]
      # collapsing shorted edge and all triangles incident to it
      # **hoping (ab) is not `boundary` because it disappears**
      # NOTE: this will remove 1 edges and 1 triangle since they will become degenerate
      if num_small==1
        shortest_edge, highest_h,(a,b,c), (ab,bc,ca)=triangle_info[:order_by_shortest_edge]
        println("Removing short edge shortest=", ab, " triangle=",[a,b,c], " other_norm=", [bc,ca])
        # VIEWTRIANGLES(V, triangles, explode=[1.8,1.8,1.8])
        return transform_triangles(V, segments, triangles, vertex_map=Dict( b=> a )) 
      end

      if num_small==0

        longest_edge, shortest_h, (a,b,c), (ab,bc,ca) = triangle_info[:order_by_longest_edge]
        if shortest_h>epsilon
          continue
        end

        # MOST COMPLICATE CASE if I want to preserve the topology
        # there is also the case that the triangle is very skeewed (its height is very small)
        # need to replace the long edge (a,b) with the shortes edge (a,c) (c,b)
        # if there is an adiancent triangle using the (a,b vertex), this triangle will be split into two triangles or collapsed depending on the configuration

        # about boundary...
        #   if [a,b] is     in segments it will replaced by two new edges
        #   if [a,b] is NOT in segment map nothing will happen 
        segment_map =Dict{Cell,Cells}([a,b] => [[a,c],[c,b]]) # ab  will disappear
        triangle_map=Dict{Cell,Cells}([a,b,c] =>[])           # abc will disappear

        for t2 in triangles
          if a in t2 && b in t2 && !(c in t2)
            adj=t2
            d=[K for K in t2 if K!=a && K!=b][1]
            @assert(d!=c)
            cd=LinearAlgebra.norm(V[:,d]-V[:,c])
            if cd<2*epsilon
              println("Removing skewed with near skewed and [a,b,c,d]=", [a,b,c,d], " , d->c since they are near each othger")
              # abd will disappear too since d almost equal to `c``
              triangle_map[[a,b,d]]=[] 
              # VIEWTRIANGLES(V, triangles, explode=[1.8,1.8,1.8])
              return transform_triangles(V, segments, triangles, segment_map=segment_map, triangle_map=triangle_map, vertex_map=Dict(d=>c))
            else
              println("Removing skewed with near big splitting the adj [a,b,c,d]=", [a,b,c,d])
              # abd will split into 2 triangles (i.e. cannot disappear since `d`` is different from `c`)
              triangle_map[[a,b,d]]=[[a,c,d],[b,c,d]] 
              # VIEWTRIANGLES(V, triangles, explode=[1.8,1.8,1.8])
              return transform_triangles(V, segments, triangles, segment_map=segment_map, triangle_map=triangle_map)
            end
          end
        end

        # no adjacency
        println("Removing skewed with no adj [a,b,c]=", [a,b,c])
        # VIEWTRIANGLES(V, triangles, explode=[1.8,1.8,1.8])
        return transform_triangles(V, segments, triangles, segment_map=segment_map, triangle_map=triangle_map)
      end
    end

  end

  # cannot simplify...
  return nothing

end

# ///////////////////////////////////////////////////////////////////
function fragment_lar_face(points_db::PointsDB, dst::Lar, src::Lar, F1::Int; debug_mode=false)

  plane_points_db = PointsDB()
  face_segments   = Cells()
  other_segments  = Cells()

  # prepare to project
  num_faces=length(src.C[:FV])
  
  fe1=src.C[:FE][F1]
  fv1=src.C[:FV][F1] 
  world_box1=bbox_create(src.V[:,fv1])
  plane=plane_create(src.V[:,fv1])
  projector=project_points3d(src.V[:,fv1], double_check=true) # scrgiorgio: remove double check 

  if debug_mode
    VIEWEDGES(src.V, [src.C[:EV][E1] for E1 in fe1])
  end

  # project face, find it bounding box
  begin
    face_proj_points=Vector{Vector{Float64}}()
    for E1 in fe1
      a,b=src.C[:EV][E1]
      proj1,proj2=[Vector{Float64}(it) for it in eachcol(projector(src.V[:,[a,b]]))]
      A=add_point(plane_points_db, proj1) 
      B=add_point(plane_points_db, proj2)
      if A!=B 
        append!(face_proj_points,[proj1, proj2])
        push!(face_segments,[A,B]) 
      end
    end
    face_proj_box=bbox_create(face_proj_points)
  end

  if debug_mode
    VIEWEDGES(get_points(plane_points_db), face_segments, title="arrange3d / frament $(F1) / projected edges")
  end

  # project other faces
  for F2 in 1:num_faces
    
    if (F1==F2) 
      continue 
    end

    # quick discard not intersecting faces
    fv2, fe2=src.C[:FV][F2], src.C[:FE][F2]
    world_box2=bbox_create(src.V[:,fv2])
    if !(bbox_intersect(world_box1, world_box2))
      continue
    end

    # find intersection on the main face
    begin
      proj_intersections=Vector{Vector{Float64}}()
      for E2 in fe2
        a,b = src.C[:EV][E2]
        p1,p2=[Vector{Float64}(it) for it in eachcol(src.V[:,[a,b]])]
        hit, t = plane_ray_intersection(p1, normalized(p2-p1), plane)
        if !isnothing(hit)
          @assert(t>0)
          # must cross the plane (be on opposite sides of the plane)
          if t<LinearAlgebra.norm(p2-p1)
            hit_proj=projector(hcat(hit))[:,1]
            push!(proj_intersections, hit_proj)
          end
        end
      end
    end

    # each face can intersect in 2 points ???! not sure for any face here
    # am i forcing convexy here?
    # @assert(length(proj_intersections)==0 || length(proj_intersections)==2)
    # println(F1, " ", F2," ", proj_intersections)
    if length(proj_intersections)==2
      a,b=proj_intersections
      other_proj_box=bbox_create([a,b])
      if bbox_intersect(face_proj_box, other_proj_box)
        A=add_point(plane_points_db, a)
        B=add_point(plane_points_db, b)
        if A!=B 
          push!(other_segments, [A,B]) 
        end
      end
    end

  end

  plane_points        = get_points(plane_points_db)
  plane_points_row    = BYROW(plane_points)
  face_segments       = remove_duplicates(face_segments)
  all_segments        = remove_duplicates([face_segments; other_segments])

  tin = Triangulate.TriangulateIO()
  tin.pointlist = plane_points 
  tin.segmentlist = hcat(all_segments...) # constrained triangulation

  if debug_mode
    VIEWEDGES(tin.pointlist, tin.segmentlist, title="arrange2d_v2 / input")
  end

  # https://github.com/JuliaGeometry/Triangulate.jl/blob/8e0fc2c0fb58ffb3ae351afcca5da375522d2e84/docs/src/triangle-h.md?plain=1#L193
  # -p Triangulates a Planar Straight Line Graph, 
  # -Q for quiet
  # -X  No exact arithmetic.
  # -q  Quality mesh generation by Delaunay refinement 
  # -D  Conforming Delaunay triangulatio
  # -S  Specifies the maximum number of Steiner points
  # -Y  No new vertices on the boundary.  This switch is useful when the mesh boundary must be preserved so that it conforms to some adjacent mesh
  #     USe `-YY' to prevent all segment splitting, including internalboundaries.
  tout=nothing
  while true
    try
      tout, __ = Triangulate.triangulate("pQ", tin) 
      break
    catch TriangulateError
      println("WARNING Triangulate failed, perturing points")
      tin.pointlist=ToPoints([p + rand(2) * LAR_ARRANGE2D_PERTURBATION for p in eachcol(plane_points)])
    end
  end

  if debug_mode
    VIEWTRIANGLES(tout.pointlist, tout.trianglelist, explode=[1.0,1.0,1.0], title="arrange2d_v2 / output / triangles")
  end

  plane_points, Vinput=tout.pointlist, plane_points
  segments  = simplify_cells([Cell([a,b])   for (a,b  ) in eachcol(tout.segmentlist)])
  triangles = simplify_cells([Cell([a,b,c]) for (a,b,c) in eachcol(tout.trianglelist)])
  tout=nothing

  # remove small triangles (it can happen if I have segment very near to the boundary like for arrange3d function)
  begin
    while true
      #@show(segments)
      #@show(triangles)
      next=remove_small_or_skewed_triangles(plane_points, segments, triangles)
      if isnothing(next)  break  end
      plane_points, segments,triangles=next
    end
  end

  # I want to keep only vertices 
  #   - coming from user input
  #   - coming from the intersection of 2 EVs
  # I do NOT want spurios points on the boundary due to the triangulation (I am getting some of them even if I use -YY)
          
  begin
    
    Vint=Vector{Vector{Float64}}()
    for (E1,(a,b)) in enumerate(all_segments)
      p1,p2=Vinput[:,a],Vinput[:,b]
      for (E2,(c,d)) in enumerate(all_segments)
        if E1==E2 continue end
        p3,p4=Vinput[:,c],Vinput[:,d]
        inter=segment_intersection(p1,p2,p3,p4, LAR_ARRANGE2D_SMALL_TRIANGLES_ERR)
        println("Vint $(a) $(b) $(c) $(d) $(inter)")
        if !isnothing(inter)
          push!(Vint,Vector{Float64}([it for it in inter]))
        end
      end
    end
    Vint=ToPoints(Vint)

    keeps=Set()
    for (I, p) in enumerate(eachcol(plane_points))

      keep=false
      for (J,q) in enumerate(eachcol(Vinput))
        if point_distance(p,q) < LAR_ARRANGE2D_SMALL_TRIANGLES_ERR
          println("KEEP $(I) because near to user vertices $(J)")
          keep=true
          break
        end
      end

      for (J,q) in enumerate(eachcol(Vint))
        if point_distance(p,q) < LAR_ARRANGE2D_SMALL_TRIANGLES_ERR
          println("KEEP $(I) because near to cross intersection")
          keep=true
          break
        end
      end

      if keep
        push!(keeps,I)
      else
        println("SPURIOUS ON BOUNDARY ", I, p)
      end
      
    end
  end

  # NOTE: segment list does not contain internal edges, but only "important" edges
  a2d=Lar(plane_points)

  # compute EV,FE
  begin
    a2d.C[:FE]=Cells()
    a2d.C[:EV]=Cells()

    adjacent_triangles=find_adjacents_cells(triangles, 2, segments)
    groups=find_groups_of_cells(adjacent_triangles)

    for (A, triangle_ids) in enumerate(groups)

      # each group will form a face (even non-convex, holed)
      face_boundary=Cells()
      internal_points=Vector{Vector{Float64}}() # this is needed in case classify is run
      for triangle_id in triangle_ids 
        u,v,w = triangles[triangle_id]
        push!(internal_points,(a2d.plane_points[:,u] + a2d.plane_points[:,v] + a2d.plane_points[:,w])/3.0)
        for (a,b) in [ [u,v], [v,w], [w,u] ]
          a,b = normalize_cell([a,b])
          if [a,b] in segments
            push!(face_boundary, [a,b])
          end
        end
      end
      face_boundary=simplify_cells(face_boundary)
      cycles=find_vcycles(face_boundary)

      # skip if outside
      if !isnothing(classify_for_3d)
        outside = length([p for p in internal_points if classify_for_3d(p) == "p_out"])
        inside  = length(internal_points) - outside

        keep_face=false
        if outside>0 && inside>0
          # is this the right thing to do???
          println("WARNING ambiguity #inside=", inside," #outside=", outside)
          keep_face=inside > outside 
        else
          keep_face=inside>0
        end

        if !keep_face
          continue
        end
      end

      for cycle in cycles

        loop=[a for (a,b) in cycle]

        # insert into the lar complex (note: i will simplify later!)
        loop=[v_index for v_index in loop if v_index in keeps]
        fe=Cell()
        for I in 1:length(loop)
          a = loop[I]
          b = loop[I==length(loop) ? 1 : I+1] 
          push!(a2d.C[:EV], [a,b])
          push!(fe, length(a2d.C[:EV]))
        end
        @show(fe)
        push!(a2d.C[:FE], fe)
      end
    end
  end

  a2d=SIMPLIFY(a2d)
  COMPUTE(a2d,:FV)
  CHECK(a2d)

  if debug_mode
    VIEWCOMPLEX(a2d, explode=[1.0,1.0,1.0], show=["V","FV","Vtext"], title="arrange2d_v2 / output / final")
  end

  # store to the destination lar unprojecting points
  faces3d=Cell()
  begin
    unprojected=projector(a2d.V, inverse=true)

    for (F,fe) in enumerate(a2d.C[:FE])

      # just want to keep triangles which are inside the main face
      unprojected_fe=Cell()
      for E in fe
        a,b = a2d.C[:EV][E]
        r1=round_vector(unprojected[:,a], digits=LAR_ARRANGE3D_UNPROJECT_ROUND_DIGITS)
        r2=round_vector(unprojected[:,b], digits=LAR_ARRANGE3D_UNPROJECT_ROUND_DIGITS)
        A=add_point(points_db, r1)
        B=add_point(points_db, r2)
        if A!=B
          push!(dst.C[:EV],[A,B])
          push!(unprojected_fe, length(dst.C[:EV]))
        end
      end

      # must still be a face with cycles
      push!(dst.C[:FE], unprojected_fe)
      push!(faces3d,length(dst.C[:FE]))

      n1=length(find_vcycles([dst.C[:EV][E] for E in dst.C[:FE][end]]))
      #n2=length(find_vcycles([a2d.C[:EV][E] for E in dst.C[:FE][end]]))
      #@assert(n1==n2)
      
    end
  end

  # TODO: faster
  dst.V=get_points(points_db)
  if haskey(dst.C,:FV)  delete!(dst.C,:FV)  end
  COMPUTE(dst,:FV)
  
  if debug_mode 
    VIEWCOMPLEX(SELECT(dst,faces3d), show=["FV", "Vtext"],explode=[1.0, 1.0, 1.0], title="arrange3d / frament $(F1) / 3d")
  end

end

# ///////////////////////////////////////
#  TGW 3D from here
# ///////////////////////////////////////


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

# //////////////////////////////////////////////////////////////
function arrange3d_v2(lar::Lar; debug_mode=true)

  # fragment
  begin
    points_db=PointsDB() 
    fragmented=Lar(zeros(Float64,0,0), Dict{Symbol,Cells}( :EV => Cells(), :FE => Cells()))
    for (F, fv) in enumerate(lar.C[:FV])
      println("# ////////////////// fragmenting ",F,"/",length(lar.C[:FV]))
      fragment_lar_face(points_db, fragmented, lar, F, debug_mode=debug_mode)
    end
    lar=SIMPLIFY(fragmented)
  end


  # if you had a problem with an edge, you should be able to debug here
  # debug_edge(lar,21,26)

  if debug_mode
    VIEWCOMPLEX(lar, explode=[1.0,1.0,1.0], show=["V","FV","Vtext"], title="arrange3d_v2 after all faces fragmentation")
  end  

  # will be created again later
  delete!(lar.C,:FV)

  # find cycles (note: the same cycle can appear twice, so I need to filter it out)
  begin

    # need to create the new FE (since two cycles will become the same face)
    # note: ev will be the same with the same indices
    existing_edges=Dict( normalize_cell([a,b]) => E for (E,(a,b)) in enumerate(lar.C[:EV]) )
    FE=Dict{Cell,Cycle}()
    for (F,fe) in enumerate(lar.C[:FE])
      face_edges=[lar.C[:EV][E] for E in fe]
      for cycle in find_vcycles(face_edges)
        # println("Face ",F, " has cycle ", cycle, " face_edges=",face_edges)
        fe=normalize_cell([existing_edges[normalize_cell([a,b])] for (a,b) in cycle])
        if !haskey(FE,fe)
          FE[fe]=cycle
        end
      end
    end
    lar.C[:FE], cycles=collect(keys(FE)),collect(values(FE))
  end

  COMPUTE(lar,:FV)

  lar.C[:CF]=lar_find_atoms(lar.V, cycles, debug_mode=debug_mode)
  
  if debug_mode
    VIEWCOMPLEX(lar, show=["FV","atom"], explode=[1.2,1.2,1.2], title="arrange3d_v2 showing atoms")
  end

  return lar
end
export arrange3d_v2

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










# /////////////////////////////////////////////////////////////
PointsDB=Dict{Vector{Float64},Int}
export PointsDB

# ///////////////////////////////////////////////////////////////////
function round_vector(v::Vector{Float64}; digits::Int)::Vector{Float64}
  if digits==0 return v end
  return [round(value, digits=digits) for value in v]
end
export round_vector

# ///////////////////////////////////////////////////////////////////
function add_point(db::PointsDB, p::Vector{Float64})::Int

  p=[abs(it)==0.0 ? 0.0 : it for it in p]

  if !haskey(db, p)  
    db[p]=length(db)+1 
  end
  return db[p]
end
export add_point

# ///////////////////////////////////////////////////////////////////
function get_points(db::PointsDB)::Points
	v=[Vector{Float64}() for I in 1:length(db)]
	for (pos,idx) in db
		v[idx]=pos
	end
	return hcat(v...)
end
export get_points

# ///////////////////////////////////////////////////////
function find_adjacents(cells::Matrix, needed_len::Int, uncrossable::Dict)::Dict
  tot=size(cells,2)
  ret= Dict( (id => Set{Int}() for id in 1:tot) )
  for (id1, v1) in enumerate(eachcol(cells))
    for (id2, v2) in enumerate(eachcol(cells))
      if id1==id2 continue end
      intersection=collect(intersect(Set(v1),Set(v2)))
      if length(intersection)==needed_len
        intersection=collect(sort(intersection))
        if !haskey(uncrossable,intersection)
          push!(ret[id1], id2)
          push!(ret[id2], id1)
        end
      end
    end
  end
  return ret
end

# ///////////////////////////////////////////////////////
function find_groups(adjacent::Dict)
  ret=[]
  to_assign=Set(1:length(keys(adjacent)))
  while !isempty(to_assign)
    seed=pop!(to_assign)
    group=[seed]
    push!(ret, group)
    stack=[seed]
    while !isempty(stack)
      cur=popfirst!(stack)
      for other in adjacent[cur]
        if other in group
          # all ok, aready assigned

        elseif other in to_assign
          delete!(to_assign,other)
          push!(group, other)
          push!(stack, other)
          
        else
          # internal error, show not happen
          @assert(false) 
        end
      end
    end
  end
  return ret
end

# /////////////////////////////////////////////////////////////////////
function print_matrix_by_col(name::String,value::Matrix)
  println(name); for (I,it) in enumerate(eachcol(value)) println(I," ",it) end
end

# /////////////////////////////////////////////////////////////////////
function show_edges(V::Points, EV::Cells; explode=[1.0,1.0,1.0])
  lar=Lar(V,Dict{Symbol,Cells}(:EV=>EV))
  # print_matrix_by_col("sorted_points", hcat(collect(sort([it for it in eachcol(V)]))))
  VIEWCOMPLEX(lar, explode=explode, show=["V","EV","Vtext"])
end

# /////////////////////////////////////////////////////////////////////
function show_edges(V::Points, segmentlist::Matrix; explode=[1.0,1.0,1.0])
  show_edges(V,[Cell(it) for it in eachcol(segmentlist)], explode=explode)
end

# /////////////////////////////////////////////////////////////////////
function show_triangles(V::Points, triangles::Matrix; explode=[1.0,1.0,1.0])
  lar=Lar(V, Dict{Symbol,Cells}(:EV => Cells(),:FE => Cells()))
  for (u,v,w) in eachcol(triangles)
    E=length(lar.C[:EV])
    append!(lar.C[:EV], [[u,v],[v,w],[w,u]])
    push!(lar.C[:FE], [E+1,E+2,E+3])
  end
  compute_FV(lar)
  VIEWCOMPLEX(lar, explode=explode, show=["V", "EV", "FV", "Vtext"])
end

# /////////////////////////////////////////////////////////////////////
function arrange2d_experimental(V::Points, EV::Cells; debug_mode=false, classify=nothing)

  tin = Triangulate.TriangulateIO()
  tin.pointlist = V 

  # constrained triangulation
  tin.segmentlist = hcat(EV...) 

  if debug_mode
    println("# tin")
    print_matrix_by_col("# pointlist", tin.pointlist )
    print_matrix_by_col("# EV"       , tin.segmentlist)
    show_edges(tin.pointlist, tin.segmentlist)
  end

  # https://github.com/JuliaGeometry/Triangulate.jl/blob/8e0fc2c0fb58ffb3ae351afcca5da375522d2e84/docs/src/triangle-h.md?plain=1#L193
  # -p Triangulates a Planar Straight Line Graph, 
  # -Q for quiet
  # -X  No exact arithmetic.
  # -q  Quality mesh generation by Delaunay refinement 
  # -D  Conforming Delaunay triangulatio
  #   -S  Specifies the maximum number of Steiner points
  tout=nothing
  while true
    try
      tout, __ = Triangulate.triangulate("Qp", tin) 
      break
    catch TriangulateError
      println("WARNING Triangulate failed to perturing points")
      tin.pointlist=ToPoints([p + LAR_FRAGMENT_ERR*rand(2) for p in eachcol(V)])
    end
  end

  if debug_mode
    println("# tout")
    print_matrix_by_col("# pointlist"   , tout.pointlist   )
    print_matrix_by_col("# segmentlist" , tout.segmentlist )
    print_matrix_by_col("# trianglelist", tout.trianglelist)

    # show_edges(tout.pointlist, tout.segmentlist, explode=[1.2,1.2,1.2])
    # show_triangles(tout.pointlist, tout.trianglelist, explode=[1.2,1.2,1.2])
  end

  # NOTE: segment list does not contain internal edges, but only "important" edges
  ret=Lar(tout.pointlist, Dict{Symbol,Cells}())

  # EV (note: segmentlist contains only boundary edges == edges that cannot be crossed; all other edges are internal)
  begin
    ret.C[:EV]=Cells()
    edges = Dict{Vector{Int},Int}()
    for (E,(a,b))  in enumerate(eachcol(tout.segmentlist))
      edges[ sort([a,b]) ]=E
      push!(ret.C[:EV], [a,b])
    end
  end

  # FE
  begin
    ret.C[:FE]=Cells()
    adj=find_adjacents(tout.trianglelist, 2, edges)
    groups=find_groups(adj)
    for (A, triangle_ids) in enumerate(groups)
      # each group will form a face (even non-convex, holed)
      fe=Cell()
      # this is needed for containment testing
      internal_points=Vector{Vector{Float64}}() 
      for triangle_id in triangle_ids
        (u,v,w) = tout.trianglelist[:,triangle_id]
        push!(internal_points,(ret.V[:,u]+ret.V[:,v]+ret.V[:,w])/3.0)
        for (a,b) in [collect(sort([u,v])),collect(sort([v,w])),collect(sort([w,u]))]
          if haskey(edges,[a,b])
            push!(fe, edges[ [a,b] ])
          end
        end
      end

      fe=remove_duplicates(fe)
      if (length(fe)>=3)
        
        # keep only internal
        if !isnothing(classify)
          outside = length([p for p in internal_points if classify(p) == "p_out"])
          inside  = length(internal_points) - outside
          if !keep_face(inside, outside)
            continue
          end
        end

        push!(ret.C[:FE], fe)

      end
    end
  end

  compute_FV(ret)

  ret=SIMPLIFY(ret)

  if debug_mode
    VIEWCOMPLEX(ret, explode=[1.2,1.2,1.2], show=["V","EV","FV","Vtext"])
  end

  return ret
end
export arrange2d_experimental

# /////////////////////////////////////////////////////////////////////
function arrange2d_experimental(lar::Lar)
  return arrange2d_experimental(lar.V,lar.C[:EV])
end

# ///////////////////////////////////////////////////////////////////
function keep_face(inside::Int, outside::Int)::Bool
  if inside==0 && outside==0
    return false

  elseif inside>0 && outside==0
    return true

  elseif outside>0 && inside==0
    return false

  else
    # is this the right thing to do???
    println("WARNING ambiguity #inside=", inside," #outside=", outside)
    return inside > outside 
  end
end

# ///////////////////////////////////////////////////////////////////
function fragment_lar_face(points_db::PointsDB, dst::Lar, src::Lar, F1::Int; debug_mode=false)

  # debug_mode=(F1==4)

  plane_points_db = PointsDB()
  face_segments   = Cells()
  other_segments  = Cells()

  # prepare to project
  num_faces=length(src.C[:FV])
  fv1, fe1=src.C[:FV][F1], src.C[:FE][F1]
  world_box1=bbox_create(src.V[:,fv1])
  plane=plane_create(src.V[:,fv1])
  projector=project_points3d(src.V[:,fv1], double_check=true) # scrgiorgio: remove double check 

  #if debug_mode
  #  println("Fragmenting face ",F1)
  #  print_matrix_by_col("V", src.V[:,fv1])
  #  println("plane",plane)
  #  show_edges(src.V[:,fv1], [src.C[:EV][E1] for E1 in fe1])
  #end

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
    show_edges(get_points(plane_points_db), face_segments)
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

            # I need to snap to existing vertices,to avoid to insert too many vertices
            # I am doing this only for "other" faces
            if true
              nearest=nothing
              for (existing,index) in plane_points_db
                d=LinearAlgebra.norm(hit_proj-existing)
                if d<LAR_FRAGMENT_ERR && (isnothing(nearest) || d<nearest)
                  nearest, hit_proj = d, existing
                end
              end
            end
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

  a2d=arrange2d_experimental(plane_points, all_segments, debug_mode=debug_mode, classify = pos -> classify_point(pos, plane_points_row, face_segments) )

  if F1==4
    @show(a2d)
    VIEWCOMPLEX(a2d,show=["FV", "Ftext", "Vtext"],explode=[2.2,1.2,1.2])
  end

  # store to the destination lar unprojecting points
  begin
    unprojected=projector(a2d.V, inverse=true)
    for (F,fe) in enumerate(a2d.C[:FE])

      # just want to keep triangles which are inside the main face
      unprojected_fe=Cell()
      fv=Cell()
      for E in fe
        a,b = a2d.C[:EV][E]
        r1=round_vector(unprojected[:,a], digits=LAR_FRAGMENT_DIGITS)
        r2=round_vector(unprojected[:,b], digits=LAR_FRAGMENT_DIGITS)
        A=add_point(points_db, r1)
        B=add_point(points_db, r2)
        if A!=B
          push!(dst.C[:EV],[A,B])
          push!(unprojected_fe,length(dst.C[:EV]))
          push!(fv,A)
          push!(fv,B)
        end
      end
      fv=remove_duplicates(fv)
      println(F1, " " ,fv)

      # must still be a face
      if length(unprojected_fe)>=3
        push!(dst.C[:FE],unprojected_fe)
        #if remove_duplicates(fv)==remove_duplicates([1, 8, 25, 27, 28])
        #  @show(F1, " " ,length(dst.C[:FE]), " ",fv)
        #  VIEWCOMPLEX(dst, show=["FV", "Ftext", "Vtext"], explode=[1.0,1.0,1.0], face_color=TRANSPARENT)
        #  aaa()
       # end

      end
    end
  end

end
export fragment_lar_face



# ///////////////////////////////////////////////////////////
function fragment_lar(lar::Lar; debug_mode=false)
  points_db=PointsDB() 
  ret=Lar(zeros(Float64,0,0), Dict{Symbol,Cells}( :EV => Cells(), :FE => Cells()))
  num_faces=length(lar.C[:FV])
  for (F, fv) in enumerate(lar.C[:FV])
    println("fragmenting face ",F,"/",num_faces)
    fragment_lar_face(points_db, ret, lar, F, debug_mode=false)
  end
  println("All faces fragmented")
  ret.V=get_points(points_db)
  compute_FV(ret)
  return SIMPLIFY(ret)
end
export fragment_lar


# ///////////////////////////////////////
#  TGW 3D from here
# ///////////////////////////////////////

const Cycle=Vector{Vector{Int}}
const Cycles=Vector{Cycle}

function reverse_cycle(value::Cycle)::Cycle
  return [[b, a] for (a, b) in reverse(value)]
end



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

  #v=remove_duplicates([abs(it) for it in [56, -79, 81, 30, 41, 43, -67, -70, -77, 80, -43, 84, -78, -25, 65, -62, -83, -36, -64, -41]])
  #VIEWCOMPLEX(SELECT(lar, v), show=["FV", "Ftext"], explode=[1.2,1.2,1.2], face_color=TRANSPARENT)

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
function arrange3d_experimental(lar::Lar; debug_mode=false)

  # problema qui sulla faccia 17 torna una cosa senza senso

  lar=fragment_lar(lar, debug_mode=debug_mode)

  # topology check: each edge should be shared by >=2 faces
  begin

  end

  # if you had a problem with an edge
  if false
    #@show(lar)
    #for (I,fv) in enumerate(lar.C[:FV])
    #  println("fv ", fv, " # ",I)
    #end
    #a,b=[28,27]

    #sel=[a,b]
    #for (P,p) in enumerate(eachcol(lar.V))
    #  if LinearAlgebra.norm(p - lar.V[:,a])<0.1 || LinearAlgebra.norm(p - lar.V[:,b])<0.1
    #    push!(sel,P)
    #  end
    #end
    #@show(sel)
    #sel=remove_duplicates([F for (F,fv) in enumerate(lar.C[:FV]) if length(intersect(Set(sel),Set(fv)))>0])
    #@show(sel)
    #sel=[17,35]
    #VIEWCOMPLEX(SELECT(lar,sel), show=["FV", "Ftext", "Vtext"], explode=[1.0,1.0,1.0], face_color=TRANSPARENT)
  end

  if debug_mode
    @show(lar)
    VIEWCOMPLEX(lar, explode=[1.2,1.2,1.2], show=["V","EV","FV","Vtext", "Ftext"])
  end  

  # will be created again later
  delete!(lar.C,:FV)

  # find cycles (note: the same cycle can appear twice, so I need to filter it out)
  begin

    # need to create the new FE (since two cycles will become the same face)
    # note: ev will be the same with the same indices
    existing_edges=Dict( sort([a,b]) => E for (E,(a,b)) in enumerate(lar.C[:EV]) )
    FE=Dict{Cell,Cycle}()
    for (F,fe) in enumerate(lar.C[:FE])
      face_edges=[lar.C[:EV][E] for E in fe]
      for cycle in find_vcycles(face_edges)
        println("Face ",F, " has cycle ", cycle)

        fe=remove_duplicates([existing_edges[sort([a,b])] for (a,b) in cycle])
        if !haskey(FE,fe)
          FE[fe]=cycle
        end
      end
    end
    lar.C[:FE],cycles=collect(keys(FE)),collect(values(FE))
  end

  compute_FV(lar)


  # @show(lar)
  # VIEWCOMPLEX(lar,show=["FV","Ftext","Vtext"], explode=[1.0,1.0,1.0],face_color=TRANSPARENT)
  #piece=SELECT(lar,[41,42, 67, 80, 83])
  #VIEWCOMPLEX(piece,show=["FV","Ftext","Vtext"], explode=[1.0,1.0,1.0],face_color=TRANSPARENT)

  lar.C[:CF]=lar_find_atoms(lar.V, cycles, debug_mode=debug_mode)
  
  if debug_mode
    VIEWCOMPLEX(lar, show=["FV","atom"], explode=[1.2,1.2,1.2])
  end

  return lar
end
export arrange3d_experimental

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
      print("OK guess_boundary_faces #attempt=",attempt, " distances=", distances)
      distances=sort(distances)
      return [ distances[1][end], distances[end][end]]
    end

    println("FAILED guess_boundary_faces #attempt=", attempt, " distances=", distances)

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
function arrange3d_experimental_split(lar::Lar)::Tuple{Lar,Lar}

  atoms=[cf for cf in lar.C[:CF]]

  # connected atoms
  begin
    @show(atoms)
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
function INNERS_experimental(lar::Lar)::Lar
  outers,inners = arrange3d_experimental_split(lar)
  return inners
end
export INNERS_experimental

# ////////////////////////////////////////////////////////////////
function OUTERS_experimental(lar::Lar)::Lar
  outers,inners=arrange3d_experimental_split(lar)
  return outers
end
export OUTERS_experimental






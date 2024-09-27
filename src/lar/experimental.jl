


# /////////////////////////////////////////////////////////////
PointsDB=Dict{Vector{Float64},Int}
export PointsDB

# ///////////////////////////////////////////////////////////////////
function round_vector(v::Vector{Float64}, digits::Int)::Vector{Float64}
  if digits==0 return v end
  ret=[round(value, digits=digits) for value in p]
  ret=[ret==0.0 ? 0.0 : value for value in p] # for positive negative integer
  return ret
end

# ///////////////////////////////////////////////////////////////////
function add_point(db::PointsDB, p::Vector{Float64})::Int
  if !haskey(db, p)  
    db[p]=length(db)+1 
  end
  return db[p]
end

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
  VIEWCOMPLEX(lar, explode=explode, show=["V","EV","V_text"])
end

# /////////////////////////////////////////////////////////////////////
function show_edges(V::Points, segmentlist::Matrix; explode=[1.0,1.0,1.0])
  show_edges(V,[Cell(it) for it in eachcol(segmentlist)], explode=explode)
end

# /////////////////////////////////////////////////////////////////////
function show_triangles(V::Points, triangles::Matrix{Int}; explode=[1.0,1.0,1.0])
  lar=Lar(V, Dict{Symbol,Cells}(:EV => Cells()))
  for (u,v,w) in eachcol(triangles)
    append!(lar.C[:EV], [[u,v],[v,w],[w,u]])
  end
  VIEWCOMPLEX(lar, explode=explode, show=["V", "EV", "FV", "V_text"])
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

  # -p Triangulates a Planar Straight Line Graph, 
  # -Q for quiet
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
    # show_edges(tout.pointlist, tout.segmentlist)
    # show_triangles(tout.pointlist, tout.trianglelist)
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
    VIEWCOMPLEX(ret, explode=[1.2,1.2,1.2], show=["V","EV","FV","V_text"])
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
      A=add_point(plane_points_db, proj1) # NOTE: no round here (!)
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
            # I am doing this only for other faces
            begin
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
      other_proj_box=bbox_create(proj_intersections)
      if bbox_intersect(face_proj_box, other_proj_box)
        A=add_point(plane_points_db, proj_intersections[1])
        B=add_point(plane_points_db, proj_intersections[2])
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

  # store to the destination lar unprojecting points
  begin
    unprojected=projector(a2d.V, inverse=true)
    for (F,fe) in enumerate(a2d.C[:FE])

      # just want to keep triangles which are inside the main face


      unprojected_fe=Cell()
      for E in fe
        a,b = a2d.C[:EV][E]
        A=add_point(points_db, unprojected[:,a])
        B=add_point(points_db, unprojected[:,b])
        if A!=B
          push!(dst.C[:EV],[A,B])
          push!(unprojected_fe,length(dst.C[:EV]))
        end
      end

      # must still be a face
      if length(unprojected_fe)>=3
        push!(dst.C[:FE],unprojected_fe)
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
    fragment_lar_face(points_db, ret, lar, F, debug_mode=debug_mode)
  end
  println("All faces fragmented")
  ret.V=get_points(points_db)
  compute_FV(ret)
  return SIMPLIFY(ret)
end
export fragment_lar


# /////////////////////////////////////////////////////////////////////
"""

not working (!)

function arrange3d_experimental(lar::Lar; debug_mode=false)

  # needed step: fragment all faces so that the input is PLCs (intersecting on edges)
  # otherwise it will not work
  # https://wias-berlin.de/software/tetgen/switches.d.html
  # https://wias-berlin.de/software/tetgen/plc.html
  # The definition of PLCs requires that they must be closed under taking intersections, that is two segments only can intersect at a shared point, 
  # two facets are either 
  #   - completely disjointed or 
  #   - intersecting only at shared segments or vertices  
  lar=fragment_lar(lar, debug_mode=debug_mode)

  if debug_mode
    @show(lar)
    VIEWCOMPLEX(lar, explode=[1.2,1.2,1.2], show=["V","EV","FV","V_text", "FV_text"])
  end

  facets=Cells()
  for fe in lar.C[:FE]
    edges=[lar.C[:EV][E] for E in fe]
    for cycle in find_vcycles(edges)
      push!(facets,[a for (a,b) in cycle])
    end
  end
  facets=remove_duplicates(facets)

  tin=TetGen.RawTetGenIO{Cdouble}()
  tin.pointlist=copy(lar.V)
  TetGen.facetlist!(tin, facets)

  if debug_mode
    println("# TIN")
    print_matrix_by_col("#pointlist",tin.pointlist)
    println("#facets"); for (I,it) in facets println(I, " ", it) end
  end

  # https://wias-berlin.de/software/tetgen/1.5/doc/manual/manual.pdf
  #  see Command Line Switches
  #   -Q  Quiet:  No terminal output except errors.
  #   -p  Tetrahedralizes a piecewise linear complex (PLC).
  #   -c Retains the convex hull of the PLC.
  #   -Y preserves the input surface mesh 
  #   -V Verbose: Detailed information, more terminal output.

  # perturbate a little
  # tin.pointlist=ToPoints([p+LAR_FRAGMENT_ERR*rand(3) for p in eachcol(tin.pointlist)])

  println("Running tetrahedralize...")
  tout = tetrahedralize(tin, "cpQ") 
  println("done")

  if debug_mode
    println("# TOUT")
    print_matrix_by_col("pointlist"      ,tout.pointlist       )
    print_matrix_by_col("edgelist"       ,tout.edgelist        )
    print_matrix_by_col("trifacelist"    ,tout.trifacelist     )
    print_matrix_by_col("tetrahedronlist",tout.tetrahedronlist )
  end

  # NOTE: segment list does not contain internal edges, but only "important" edges
  ret=Lar(tout.pointlist, Dict{Symbol,Cells}(
    :EV => Cells(),
    :FE => Cells(),
    :CF => Cells()))

  # edges
  begin
    edges = Dict{Vector{Int},Int}()
    for (E,(a,b))  in enumerate(eachcol(tout.edgelist))
      a,b = sort([a,b])
      @assert(a!=b)
      edges[ [a,b] ]=E
      push!(ret.C[:EV], [a,b])
    end
    if debug_mode
      @show(ret)
      VIEWCOMPLEX(ret)
    end
  end

  # faces
  begin
    faces = Dict{Vector{Int},Int}() # from (u,v,w) to face id 
    for triangles in find_groups(find_adjacents(tout.trifacelist, 2, edges))
      push!(ret.C[:FE], Cell())
      for (u,v,w) in [tout.trifacelist[:,T] for T in triangles]
        faces[sort([u,v,w])]=length(ret.C[:FE])
        for (a,b) in [sort([u,v]),sort([v,w]),sort([w,u])]
          if haskey(edges,[a,b])
            push!(ret.C[:FE][end],edges[[a,b]])
          end
        end
      end
      @assert(length(ret.C[:FE][end])>=3)
    end
    compute_FV(ret)
    if debug_mode
      VIEWCOMPLEX(ret, explode=[1.2,1.2,1.2])
    end
  end

  # CF
  begin
    for tets in find_groups(find_adjacents(tout.tetrahedronlist, 3, faces))
      push!(ret.C[:CF],Cell())
      for (u,v,w,z) in [tout.tetrahedronlist[:,T] for T in tets]
        for (a,b,c) in [sort([u,v,w]),sort([u,v,z]),sort([u,w,z]),sort([v,w,z])]
          if haskey(faces,[a,b,c])
            push!(ret.C[:CF][end], faces[ [a, b, c] ])
          end
        end
      end
      if length(ret.C[:CF][end])<4
        pop!(ret.C[:CF])
      end
    end
    if debug_mode
      @show(ret)
      VIEWCOMPLEX(ret, explode=[1.9,1.9,1.9], show=["V","EV","FV","atom"])
    end
  end

  return SIMPLIFY(ret)

end
export arrange3d_experimental#
"""


"""

# attempt to find CF
# TOTALLY WRONG, since I cannot use the bounding-box or volumes to find minimal cycles
# but keeping the code around anyway


using Plasm
using DataStructures

# ///////////////////////////////////////////////
function get_adjacent_edges(lar::Lar, active_edges::Set{Int}, E::Int)::Vector{Tuple{Int,Int}}
  ret=Vector{Tuple{Int,Int}}()
  for v_index in lar.C[:EV][E]
    num_links=length([e_index for e_index in lar.C[:VE][v_index] if e_index in active_edges])
    for e_index in lar.C[:VE][v_index]
      if (e_index!=E) && (e_index in active_edges) 
        push!(ret,(e_index, num_links))
      end
    end
  end
  return ret
end

# ///////////////////////////////////////////////
function compute_edge_bbox(lar, E::Int)::Vector
  return bbox_create(lar.V[:,lar.C[:EV][E]])
end

# ///////////////////////////////////////////////
function compute_bbox_volume(box)::Float64
  m,M=box
  return prod([(b-a) for (a,b) in zip(m,M)])
end

# ///////////////////////////////////////////////
function compute_box_union(b1,b2)
  m1,M1=b1
  m2,M2=b2
  return minimum([m1,m2]),maximum([M1,M2])
end

# ///////////////////////////////////////////////
function find_boundary_edge(lar::Lar, active_edges::Set{Int}; max_attempts=100)::Int
  
  points=[]
  for E in active_edges
    fv=lar.C[:EV][E]
    append!(points,[p for p in eachcol(lar.V[:, fv])])
  end

  b1,b2=bbox_create(ToPoints(points))
  move_out = LinearAlgebra.norm(b2 - b1)

  for attempt in 1:max_attempts

    inside_bbox    = [random(b1[I],b2[I]) for I in 1:length(b1) ]
    external_point = inside_bbox + move_out * random_dir()
    ray_dir        = normalized(inside_bbox-external_point)

    num_hits=0
    first_hit=nothing
    for E in active_edges
      hit, distance=ray_face_intersection(external_point, ray_dir, lar, F)
      if !isnothing(hit) 
        num_hits+=1
        if isnothing(first_hit) || distance<first_hit[2]
          first_hit= E, distance
        end
      end
    end

    # should go from outside to outside
    if (num_hits >= 2) && (num_hits % 2) == 0
      return first_hit[1]
    end

  end

  # cannot find anyone
  @assert(false)
    
end

# ///////////////////////////////////////////////
function find_boundary_edges(result::Vector(), lar::Lar, active_edges::Set{Int}, First::Int)

  @assert(First in active_edges)

  # already added
  if First in active_edges
    return 
  end

  push!(result, First)

  for (Adj, num_links) in get_adjacent_edges(lar, active_edges, First)
    if !(Adj in active_edges) continue
    if num_links==2 # it must be boundary
      find_boundary_edges(result, lar, active_edges, Adj)
    end
  end

end


function find_boundary_edges(lar::Lar, active_edges::Set{Int}, First::Int)::Vector
  result=Vector()
  find_boundary_edges(result, lar, active_edges, E)
  return result
end


# ///////////////////////////////////////////////
function find_atoms(lar::Lar)

  lar.C[:VE]=compute_VE(lar)

  active_edges=Set(1::length(lar.C[:FE]))

  atoms=Vector{Vector{Int}}()

  while length(active_edges)>0

    First=find_boundary_edge(lar, active_edges)
    @assert(First in active_edges)
  
    visited=Set()
    pq=PriorityQueue{Float64, Vector}()
  
    box=compute_edge_bbox(lar, First)
    volume=compute_bbox_volume(box)
    push!(pq, volume => [box,First])
  
    atom=nothing
    while !isempty(pq)
      box, E= pq[peek(pq)]
      dequeue!(pq)  
      last=active_edges[end]
  
      if last in visited
        if last==First
          atom=active_edges
          break
        end
        continue
      end
  
      push!(visited, E)
  
      for (Adj, num_links) in get_adjacent_edges(lar, active_edges, E)
        if !(it in active_edges) continue
        sub_box=compute_box_union(box, compute_edge_bbox(lar, Adj)) 
        push!(pq,compute_bbox_volume(sub.box) => [sub_box, [active_edges; Adj]])
      end
  
    end
  
    @assert(!isnothing(atom))
    push!(atoms,atom)
  
    for it in find_boundary_edges(lar, active_edges, First)
      delete!(active_edges, it)
    end
  end

  delete!(lar.C[:VE])
  return atoms

end


lar=Lar()
find_atoms(lar)
"""
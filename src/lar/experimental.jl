using TetGen


# /////////////////////////////////////////////////////////////
mutable struct PointsDB
  digits::Int
	points::Dict{Vector{Float64},Int}

  # constructor
  function PointsDB(digits::Int)
    new(digits,Dict{Vector{Float64},Int}())
  end

end
export PointsDB

# ///////////////////////////////////////////////////////////////////
function add_point(db::PointsDB, p::Vector{Float64})::Int

  if db.digits!=0
	  p=[round(value, digits=db.digits) for value in p]
    p=[p==0.0 ? 0.0 : value for value in p] # for positive negative integer
  end

  if !haskey(db.points, p)  
    db.points[p]=length(db.points)+1 
  end

  idx=db.points[p]
  return idx
end


# ///////////////////////////////////////////////////////////////////
function get_points(db::PointsDB)::Points
	v=[Vector{Float64}() for I in 1:length(db.points)]
	for (pos,idx) in db.points
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
function arrange2d(V::Points,EV::Cells; keep_atom=nothing, debug_mode=false)

  triin = Triangulate.TriangulateIO()
  triin.pointlist = V 

  # constrained triangulation
  triin.segmentlist = hcat(EV...) 

  # Triangulates a Planar Straight Line Graph, Q for quiet
  triout, __ = Triangulate.triangulate("pQ", triin) 

  if debug_mode
    println("V" ); for (I,it) in enumerate(eachcol(V)) println(I, " ",it) end
    println("EV"); for (I,it) in enumerate(EV)         println(I, " ",it) end
    println("pointlist"   );for (I,it) in enumerate(eachcol(triout.pointlist   )) println("  ",I," ",it) end
    println("segmentlist" );for (I,it) in enumerate(eachcol(triout.segmentlist )) println("  ",I," ",it) end
    println("trianglelist");for (I,it) in enumerate(eachcol(triout.trianglelist)) println("  ",I," ",it) end
  end

  # NOTE: segment list does not contain internal edges, but only "important" edges
  ret=Lar(triout.pointlist, Dict{Symbol,Cells}())

  # show triangles
  if debug_mode
    tmp=Lar(triout.pointlist, Dict{Symbol,Cells}(:EV => Cells()))
    for (u,v,w) in eachcol(triout.trianglelist)
      append!(tmp.C[:EV], [[u,v],[v,w],[w,u]])
    end
    VIEWCOMPLEX(SIMPLIFY(tmp), explode=[1.2,1.2,1.2])
  end

  # EV (note: segmentlist contains only boundary edges == edges that cannot be crossed; all other edges are internal)
  begin
    ret.C[:EV]=Cells()
    edges = Dict{Vector{Int},Int}()
    for (E,(a,b))  in enumerate(eachcol(triout.segmentlist))
      edges[ sort([a,b]) ]=E
      push!(ret.C[:EV], [a,b])
    end
    if debug_mode
      VIEWCOMPLEX(ret, explode=[1.2,1.2,1.2])
    end
  end

  # FE
  begin
    ret.C[:FE]=Cells()
    for (A, triangle_ids) in enumerate(find_groups(find_adjacents(triout.trianglelist, 2, edges)))
      fe=Cell()
      internal_points=[]
      for triangle_id in triangle_ids
        (u,v,w) = triout.trianglelist[:,triangle_id]
        push!(internal_points,(ret.V[:,u]+ret.V[:,v]+ret.V[:,w])/3.0)
        for (a,b) in [collect(sort([u,v])),collect(sort([v,w])),collect(sort([w,u]))]
          if haskey(edges,[a,b])
            push!(fe, edges[ [a,b] ])
          end
        end
      end
      if isnothing(keep_atom) || keep_atom(internal_points)
        push!(ret.C[:FE], fe)
      end
    end
  end

  # FV
  compute_FV(ret)

  return SIMPLIFY(ret)
end
export arrange2d

# /////////////////////////////////////////////////////////////////////
function arrange2d(lar::Lar)
  return arrange2d(lar.V,lar.C[:EV])
end

# ///////////////////////////////////////////////////////////////////
function fragment_face(points_db::PointsDB, dst::Lar, src::Lar, F1::Int; debug_mode=false)

  plane_points_db = PointsDB(0) # here I do not want any approximation
  face_segments   = Cells()
  other_segments  = Cells()

  num_faces=length(src.C[:FV])

  fv1=src.C[:FV][F1]
  fe1=src.C[:FE][F1]

  plane=plane_create(src.V[:,fv1])
  projector=project_points3d(src.V[:,fv1], double_check=true) # scrgiorgio: remove double check 

  if debug_mode
    @show(BYROW(src.V[:,fv1]))
    @show(plane)
  end

  for E1 in fe1
    a,b=src.C[:EV][E1]
    p1::Vector{Float64},p2::Vector{Float64} = [it for it in eachcol(src.V[:,[a,b]])]

    A=add_point(plane_points_db, projector(hcat(p1))[:,1])
    B=add_point(plane_points_db, projector(hcat(p2))[:,1])
    if A!=B push!(face_segments,[A,B]) end
  end
  
  for F2 in 1:num_faces
    
    if (F1==F2) continue end

    fv2=src.C[:FV][F2]
    fe2=src.C[:FE][F2]

    # quick discard not intersecting faces
    world_box1=bbox_create_from_points(src.V[:,fv1])
    world_box2=bbox_create_from_points(src.V[:,fv2])
    if !(bbox_intersect(world_box1, world_box2))
      continue
    end

    intersections=[]
    for E2 in fe2
      a,b=src.C[:EV][E2]
      p1::Vector{Float64},p2::Vector{Float64}=[it for it in eachcol(src.V[:,[a,b]])]
      hit, t = plane_ray_intersection(p1, normalized(p2-p1), plane)
      if !isnothing(hit)
        @assert(t>0)
        # must cross the plane
        if t<LinearAlgebra.norm(p2-p1)
          push!(intersections, hit)
        end
      end
    end

    # each face can intersect in 2 points
    # am i forcing convexy here?
    # @assert(length(intersections)==0 || length(intersections)==2)
    # println(F1, " ", F2," ", intersections)
    if length(intersections)==2
      A=add_point(plane_points_db, projector(hcat(intersections[1]))[:,1])
      B=add_point(plane_points_db, projector(hcat(intersections[2]))[:,1])
      if A!=B push!(other_segments,[A,B]) end
    end

  end

  plane_points        = get_points(plane_points_db)
  plane_points_row    = BYROW(plane_points)
  face_segments       = remove_duplicates(face_segments)
  all_segments        = remove_duplicates([face_segments; other_segments])

  # if you want to debug prjected points
  if debug_mode
    VIEWCOMPLEX(Lar(plane_points, Dict{Symbol,Cells}( :EV => all_segments)))
  end

  # the logic here is that I want to have only the fragment face, not all "around"
  function is_inside_face(centroids::Vector)
    inside,outside=0,0
    for centroid in centroids
      if classify_point(centroid, plane_points_row, face_segments) == "p_out"
        outside+=1
      else
        inside+=1
      end
    end

    if inside==0
      return false

    elseif outside==0
      return true
    else
      # both positive, cannot decide really...si this the right thing to do???
      println("   # centroids=", length(centroids), " #inside=", inside," #outside=", outside)
      return inside > outside 
    end

  end

  a2d=arrange2d(plane_points, all_segments,keep_atom=is_inside_face, debug_mode=debug_mode)

  if debug_mode
    VIEWCOMPLEX(a2d)
  end

  # store to the destination lar unprojecting points
  unprojected=projector(a2d.V, inverse=true)
  
  for fe in a2d.C[:FE]
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
    if length(unprojected_fe)>=3
      push!(dst.C[:FE],unprojected_fe)
    end
  end

end
export fragment_face



# ///////////////////////////////////////////////////////////
function fragment(lar::Lar; fragment_digits::Int=LAR_FRAGMENT_DIGITS, debug_mode=false)
  points_db=PointsDB(fragment_digits) 
  ret=Lar(zeros(Float64,0,0), Dict{Symbol,Cells}( :EV => Cells(), :FE => Cells()))
  num_faces=length(lar.C[:FV])
  for (F, fv) in enumerate(lar.C[:FV])
    println("fragmenting face ",F,"/",num_faces)
    fragment_face(points_db, ret, lar, F, debug_mode=false)
  end
  println("All faces fragmented")
  ret.V=get_points(points_db)
  compute_FV(ret)
  return SIMPLIFY(ret)
end

const LAR_FRAGMENT_DIGITS=5
export LAR_FRAGMENT_DIGITS

# /////////////////////////////////////////////////////////////////////
function arrange3d(lar::Lar; debug_mode=false, fragment_digits=LAR_FRAGMENT_DIGITS)

  # needed step: fragment all faces so that the input is PLCs (intersecting on edges)
  # otherwise it will not work
  # https://wias-berlin.de/software/tetgen/switches.d.html
  # https://wias-berlin.de/software/tetgen/plc.html
  # The definition of PLCs requires that they must be closed under taking intersections, that is two segments only can intersect at a shared point, 
  # two facets are either 
  #   - completely disjointed or 
  #   - intersecting only at shared segments or vertices  
  lar=fragment(lar, fragment_digits=fragment_digits, debug_mode=false)

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
  # https://wias-berlin.de/software/tetgen/1.5/doc/manual/manual.pdf
  #  see Command Line Switches
  #   -Q  Quiet:  No terminal output except errors.
  #   -p  Tetrahedralizes a piecewise linear complex (PLC).
  #   -c Retains the convex hull of the PLC.
  #   -Y preserves the input surface mesh 
  #   -V Verbose: Detailed information, more terminal output.
  tin=TetGen.RawTetGenIO{Cdouble}()
  tin.pointlist=copy(lar.V)
  TetGen.facetlist!(tin, facets)
  tout = tetrahedralize(tin, "cpQ") 

  if debug_mode
    println("# TIN")
    println("pointlist"      );for (I,it) in enumerate(eachcol(tin.pointlist       )) println("  ",I," ",it) end
    println("facets"         );for (I,it) in enumerate(facets                      ) println("  ",I," ",it) end
    println("")
    println("# TOUT")
    println("pointlist"      );for (I,it) in enumerate(eachcol(tout.pointlist       )) println("  ",I," ",it) end
    println("edgelist"       );for (I,it) in enumerate(eachcol(tout.edgelist        )) println("  ",I," ",it) end
    println("trifacelist"    );for (I,it) in enumerate(eachcol(tout.trifacelist     )) println("  ",I," ",it) end
    println("tetrahedronlist");for (I,it) in enumerate(eachcol(tout.tetrahedronlist )) println("  ",I," ",it) end
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
export arrange3d
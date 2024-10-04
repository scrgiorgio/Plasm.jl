

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

# /////////////////////////////////////////////////////////////////////
function arrange2d_v2(V::Points, EV::Cells; debug_mode=false, classify=nothing)

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


  # remove small triangles (it can happen if I have segment very near to the boundary like for arrange3d function)
  remove_small_triangle=LAR_FRAGMENT_ERR
  if !isnothing(remove_small_triangle) && remove_small_triangle>0

    triangles=Vector{Any}([it for it in eachcol(tout.trianglelist)])
    
    boundary_segments=Set{Vector{Int}}([collect(sort([a,b])) for (a,b) in eachcol(tout.segmentlist)])
    function is_boundary(a,b)  return collect(sort([a,b])) in boundary_segments end
    function add_boundary(a,b) push!(boundary_segments,collect(sort([a,b])))    end
    function rm_boundary(a,b)  delete!(boundary_segments,collect(sort([a,b])))  end

    T=1
    while T<=length(triangles)

      if isnothing(triangles[T])
        T+=1
        continue
      end

      a,b,c=triangles[T]

      ab=LinearAlgebra.norm(tout.pointlist[:,b]-tout.pointlist[:,a])
      bc=LinearAlgebra.norm(tout.pointlist[:,c]-tout.pointlist[:,b])
      ca=LinearAlgebra.norm(tout.pointlist[:,a]-tout.pointlist[:,c])

      # order so (a,b) is the longest edge
      __longest,(ab,bc,ca),(a,b,c)=sort([
        (ab,(ab,bc,ca),(a,b,c)),
        (bc,(bc,ca,ab),(b,c,a)),
        (ca,(ca,ab,bc),(c,a,b)),
      ])[1]

      # large enough not split and go to the next one
      if (bc+ca-ab)>LAR_FRAGMENT_ERR
        T+=1
        continue
      end

      # remove this triangle and split the other one (if there is one)
      println("Removing small triangle ", bc+ca-ab)
      triangles[T]=nothing

      # if (a,b) is on the boundary, (c,a) (b,c) will become boundary
      if is_boundary(a,b)
        rm_boundary(a,b)
        add_boundary(a,c)
        add_boundary(b,c)
      end

      # fix the triangle sharing the edge (a,b)
      # split adjacent triangle into 2 triangles 
      # NOTE: new edge (c,d) is not a boundary segment (i.e. i can move from one triangle to the other)
      Adj=[Adj for Adj in 1:length(triangles) if (Adj != T) && !isnothing(triangles[Adj]) && (a in triangles[Adj]) && (b in triangles[Adj])]
      @assert(length(Adj)<=1)
      if length(Adj)!=0
        Adj=Adj[1]
        d=[v_index for v_index in triangles[Adj] if v_index!=a && v_index!=b]
        @assert(length(d)==1)
        d=d[1]
        triangles[Adj]=nothing
        push!(triangles, [a,c,d])
        push!(triangles, [b,c,d])
      end

      T+=1

    end 

    # @show(triangles)
    triangles=[Vector{Int}(it) for it in triangles if !isnothing(it)]
    tout.trianglelist = hcat(triangles...) 
    tout.segmentlist  = hcat(collect(boundary_segments)...) 

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

  # sanity check
  begin

    # each edge should be have 1 or 2 faces
    begin
      get_faces=Dict{Int,Set{Int}}()
      for (F,fe) in enumerate(ret.C[:FE])
        for E in fe
          if !haskey(get_faces,E) get_faces[E]=Set{Int}() end
          push!(get_faces[E],F)
        end
      end
      @assert(all([length(it) in [1,2] for it in values(get_faces)]))
    end

    # one face should have at least 3 edges
    begin
      for fe in ret.C[:FE]
        @assert(length(fe)>=3)
      end
    end

  end

  if debug_mode
    VIEWCOMPLEX(ret, explode=[1.2,1.2,1.2], show=["V","EV","FV","Vtext"])
  end

  return ret
end

# /////////////////////////////////////////////////////////////////////
function arrange2d_v2(lar::Lar)
  return arrange2d_v2(lar.V,lar.C[:EV])
end
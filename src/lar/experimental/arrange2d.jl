

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


  # remove triangles too small (it can happen if I have segment very near to the boundary)
  remove_small_triangle=LAR_FRAGMENT_ERR
  if !isnothing(remove_small_triangle) && remove_small_triangle>0

    triangles=[it for it in eachcol(triangles)]
    
    boundary_segments=Set([sorted([a,b]) for (a,b) in eachcol(segmentlist)])
    function is_boundary(a,b) return haskey(boundary_segments,sort([a,b])) end
    function add_boundary(a,b) add!(boundary_segments,sort([a,b])) end
    function rm_boundary(a,b) delete(!boundary_segments,sort([a,b])) end


    T=1
    while T<=length(triangles)
      
      tpoints= [it for it in eachcol(tout.pointlist[:,triangles[T])]
      (da,db,dc),(a,b,c)=sort([
        (LinearAlgebra.norm(p[2]-p[1]),(tpoints[1],tpoints[2],tpoints[3])),
        (LinearAlgebra.norm(p[3]-p[2]),(tpoints[2],tpoints[3],tpoints[1])),
        (LinearAlgebra.norm(p[1]-p[3]),(tpoints[3],tpoints[1],tpoints[2])),
      ])[1]

      # large enough not split and go to the next one
      if (db+dc-da)>LAR_FRAGMENT_ERR
        T+=1
        continue
      end

      # remove this triangle and split the other one (if there is one)
      println("Removing small triangle ", db+dc-da)
      triangles[T]=nothing

      # if (a,b) is on the boundary, (c,a) (b,c) will become boundary
      if is_boundary(a,b):
        rm_boundary(a,b)
        add_boundary(a,c)
        add_boundary(b,c)
      end

      # fix the triangle sharing the edge (a,b)
      for Other in 1:length(triangles)
        if S!=T
          other=triangles[Other]
          if a in other && b in other
            d=[it for it in other if it!=a && it!=b];assert(length(d)==1)
            d=d[1]
            # split it into 2 triangles 
            # NOTE: new edge (c,d) is not a boundary segment (i.e. i can move from one triangle to the other)
            triangles[Other]=nothing
            push!(triangles, [a,c,d])
            push!(triangles, [b,c,d])
          end
        end
      end

    end 

    tout.trianglelist = hcat([it for it in triangles is !isnothing(it)]...) 
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
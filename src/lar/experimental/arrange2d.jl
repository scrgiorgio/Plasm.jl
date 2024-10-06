

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
function find_adjacents(triangles::Cells, uncrossable::Set)::Dict{Int,Set{Int}}
  tot=length(triangles)
  ret= Dict( (id => Set{Int}() for id in 1:tot) )
  for (id1, v1) in enumerate(triangles)
    for (id2, v2) in enumerate(triangles)
      if id1==id2 
        continue 
      end
      intersection=Cell(collect(intersect(Set(v1),Set(v2))))
      if length(intersection)==2
        intersection=normalize_cell(intersection)
        if !(intersection in uncrossable)
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
function show_edges(V::Points, EV::Cells; explode=[1.0,1.0,1.0], title="LAR")
  lar=Lar(V,Dict{Symbol,Cells}(:EV=>EV))
  # print_matrix_by_col("sorted_points", hcat(collect(sort([it for it in eachcol(V)]))))
  VIEWCOMPLEX(lar, explode=explode, show=["V","EV","Vtext"], title=title)
end

# /////////////////////////////////////////////////////////////////////
function show_edges(V::Points, segmentlist::Matrix; explode=[1.0,1.0,1.0], title="LAR")
  show_edges(V,[Cell(it) for it in eachcol(segmentlist)], explode=explode, title=title)
end

# /////////////////////////////////////////////////////////////////////
function show_triangles(V::Points, triangles::Cells; explode=[1.0,1.0,1.0], title="LAR")
  lar=Lar(V, Dict{Symbol,Cells}(:EV => Cells(),:FE => Cells()))
  for (u,v,w) in triangles
    E=length(lar.C[:EV])
    append!(lar.C[:EV], [[u,v],[v,w],[w,u]])
    push!(lar.C[:FE], [E+1,E+2,E+3])
  end
  compute_FV(lar)
  VIEWCOMPLEX(lar, explode=explode, show=["V", "EV", "FV", "Vtext"], title=title)
end

function show_triangles(V::Points, triangles::Matrix; explode=[1.0,1.0,1.0], title="LAR")
  return show_triangles(V,[Cell(it) for it in eachcol(triangles)], explode=explode, title=title)
end



# /////////////////////////////////////////////////////////////////////
function get_triangle_info(V::Points, a,b,c)::Dict

  ab=LinearAlgebra.norm(V[:,b]-V[:,a])
  bc=LinearAlgebra.norm(V[:,c]-V[:,b])
  ca=LinearAlgebra.norm(V[:,a]-V[:,c])

  s = (ab + bc + ca) / 2
  area=sqrt(s * (s - ca) * (s - bc) * (s - ab))

  # sort by base ASC, height DESC
  tmp=collect(sort([
    (ab, 2*area/ab, (a,b,c), (ab, bc, ca)),
    (bc, 2*area/bc, (b,c,a), (bc, ca, ab)),
    (ca, 2*area/ca, (c,a,b), (ca, ab, bc)),
  ]))

  return Dict(
    :order_by_shortest_edge => tmp[1],
    :order_by_middle_edge   => tmp[2],
    :order_by_longest_edge  => tmp[3],
    :edge_lengths => [ab, bc, ca]
  )

end

# /////////////////////////////////////////////////////////////////////
function check_not_boundary(V::Points,segments::Cells, v)
  v=Set([normalize_cell(it) for it in v] )
  for segment in segments
    @assert(!(normalize_cell(segment) in v))
  end
end

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

      triangle_info=get_triangle_info(V, a,b,c)
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
      # check_not_boundary(segments, [[a,b],[b,c],[c,a]]) COULD BE ?
      if num_small==3
        println("Removing triangle because too small overall triangle=", [a,b,c], " edge_lengths=", edge_lengths)
        # show_triangles(V, triangles, explode=[1.8,1.8,1.8])
        return transform_triangles(V, segments, triangles, vertex_map=Dict( b=> a, c=>a )) 
      end

      # [OK]
      # collapsing triangle and its edges in a point
      # 3rd edge has module in the range [epsilon,2epsilon] and so I think it can be removed as well
      # **hoping (ab),(bc),(c,d) are not `boundaries` because they disappear**
      # NOTE: this will remove the 3 edges and 3 triangles since they will become degenerate
      # check_not_boundary(segments, [[a,b],[b,c],[c,a]]) COULD BE ?
      if num_small==2
        println("Removing triangle because 2 small edges and 1 slightly more than small=", [a,b,c], " edge_lengths=", edge_lengths)
        # show_triangles(V, triangles, explode=[1.8,1.8,1.8])
        return transform_triangles(V, segments, triangles, vertex_map=Dict( b=> a, c=>a )) 
      end

      # [OK]
      # collapsing shorted edge and all triangles incident to it
      # **hoping (ab) is not `boundary` because it disappears**
      # NOTE: this will remove 1 edges and 1 triangle since they will become degenerate
      # check_not_boundary(segments, [[a,b]]) COULD BE?
      if num_small==1
        shortest_edge, highest_h,(a,b,c), (ab,bc,ca)=triangle_info[:order_by_shortest_edge]
        println("Removing short edge shortest=", ab, " triangle=",[a,b,c], " other_norm=", [bc,ca])
        # show_triangles(V, triangles, explode=[1.8,1.8,1.8])
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
              # show_triangles(V, triangles, explode=[1.8,1.8,1.8])
              return transform_triangles(V, segments, triangles, segment_map=segment_map, triangle_map=triangle_map, vertex_map=Dict(d=>c))
            else
              println("Removing skewed with near big splitting the adj [a,b,c,d]=", [a,b,c,d])
              # abd will split into 2 triangles (i.e. cannot disappear since `d`` is different from `c`)
              triangle_map[[a,b,d]]=[[a,c,d],[b,c,d]] 
              # show_triangles(V, triangles, explode=[1.8,1.8,1.8])
              return transform_triangles(V, segments, triangles, segment_map=segment_map, triangle_map=triangle_map)
            end
          end
        end

        # no adjacency
        println("Removing skewed with no adj [a,b,c]=", [a,b,c])
        # show_triangles(V, triangles, explode=[1.8,1.8,1.8])
        return transform_triangles(V, segments, triangles, segment_map=segment_map, triangle_map=triangle_map)
      end
    end

  end

  # cannot simplify...
  return nothing

end

# /////////////////////////////////////////////////////////////////////
function point_distance(a, b)::Float64
  return LinearAlgebra.norm(a-b)
end

# /////////////////////////////////////////////////////////////////////
# function point_is_between(a, b, c, epsilon):Bool
#   return abs(point_distance(a,c) + point_distance(c,b) - point_distance(a,b)) < epsilon
# end

# /////////////////////////////////////////////////////////////////////
#  see https://gist.github.com/kylemcdonald/6132fc1c29fd3767691442ba4bc84018
function segment_intersection(p1, p2, p3, p4, epsilon::Float64)

  (x1,y1),(x2,y2),(x3,y3),(x4,y4) = p1,p2,p3,p4
  
  # parallel
  denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1)
  if abs(denom) < epsilon
    return nothing
  end

  # out of range
  ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom
  if (ua < -epsilon) || (ua > 1+epsilon) return nothing end

  # out of range
  ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom
  if (ub < -epsilon) || (ub > 1+epsilon) return nothing end

  return (
    x1 + ua * (x2-x1),
    y1 + ua * (y2-y1))
end

    


# /////////////////////////////////////////////////////////////////////
function arrange2d_v2(V::Points, EV::Cells; debug_mode=false, classify_for_3d=nothing)

  tin = Triangulate.TriangulateIO()
  tin.pointlist = V 

  # constrained triangulation
  tin.segmentlist = hcat(EV...) 

  if debug_mode
    #println("# tin")
    #print_matrix_by_col("# pointlist", tin.pointlist )
    #print_matrix_by_col("# EV"       , tin.segmentlist)
    show_edges(tin.pointlist, tin.segmentlist, title="arrange2d_v2 / input")
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
      tin.pointlist=ToPoints([p + rand(2) * LAR_ARRANGE2D_TRIANGLE_PERTURBATION for p in eachcol(V)])
    end
  end

  if debug_mode
    #println("# tout")
    #print_matrix_by_col("# pointlist"   , tout.pointlist   )
    #print_matrix_by_col("# segmentlist" , tout.segmentlist )
    #print_matrix_by_col("# trianglelist", tout.trianglelist)
    # show_edges(V, segments, explode=[1.0,1.0,1.0], title="arrange2d_v2 / output / segments")
    show_triangles(tout.pointlist, tout.trianglelist, explode=[1.0,1.0,1.0], title="arrange2d_v2 / output / triangles")
  end

  V, Vinput=tout.pointlist, V
  segments::Cells  =[Cell([a,b])   for (a,b  ) in eachcol(tout.segmentlist)]
  triangles::Cells =[Cell([a,b,c]) for (a,b,c) in eachcol(tout.trianglelist)]
  tout=nothing

  # remove small triangles (it can happen if I have segment very near to the boundary like for arrange3d function)
  begin
    while true
      #@show(segments)
      #@show(triangles)
      next=remove_small_or_skewed_triangles(V, segments, triangles)
      if isnothing(next)  break  end
      V, segments,triangles=next
    end
  end

  # I want to keep only vertices 
  #   - coming from user input
  #   - coming from the intersection of 2 EVs
  # I do NOT want spurios points on the boundary due to the triangulation (I am getting some of them even if I use -YY)
          
  begin
    
    # 
    Vint=Vector{Vector{Float64}}()
    for (E1,(a,b)) in enumerate(EV)
      p1,p2=Vinput[:,a],Vinput[:,b]
      for (E2,(c,d)) in enumerate(EV)
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
    for (I, p) in enumerate(eachcol(V))

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
  ret=Lar(V)

  # compute EV,FE
  begin
    ret.C[:FE]=Cells()
    ret.C[:EV]=Cells()

    all_boundaries = Set{Cell}([ normalize_cell([a,b]) for (a,b) in segments ])
    adjacent_triangles::Dict{Int,Set{Int}}=find_adjacents(triangles, all_boundaries)
    groups=find_groups(adjacent_triangles)

    for (A, triangle_ids) in enumerate(groups)

      # each group will form a face (even non-convex, holed)
      face_boundary=Cells()
      internal_points=Vector{Vector{Float64}}() # this is needed in case classify is run
      for triangle_id in triangle_ids 
        u,v,w = triangles[triangle_id]
        push!(internal_points,(ret.V[:,u] + ret.V[:,v] + ret.V[:,w])/3.0)
        for (a,b) in [ [u,v], [v,w], [w,u] ]
          a,b = normalize_cell([a,b])
          if [a,b] in all_boundaries
            push!(face_boundary, [a,b])
          end
        end
      end
      face_boundary=simplify_cells(face_boundary)
      cycles=find_vcycles(face_boundary)
      @show(face_boundary)
      @show(cycles)

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

        loop=[a for (a,b) in cycles[1]]
        println("Loop ",loop)

        # insert into the lar complex (note: i will simplify later!)
        loop=[v_index for v_index in loop if v_index in keeps]
        fe=Cell()
        for I in 1:length(loop)
          a = loop[I]
          b = loop[I==length(loop) ? 1 : I+1] 
          push!(ret.C[:EV], [a,b])
          push!(fe, length(ret.C[:EV]))
        end
        @show(fe)
        push!(ret.C[:FE], fe)
      end
    end
  end

  ret=SIMPLIFY(ret)

  # compute FV
  compute_FV(ret)

  # sanity check
  begin

    # eatch edge should have 2 vertives
    for ev in ret.C[:EV]
      @assert(length(ev)==2)
    end

    # each edge should be have 1 or 2 faces
    begin
      count_faces=Dict{Int,Set{Int}}()
      for (F,fe) in enumerate(ret.C[:FE])
        for E in fe
          if !haskey(count_faces,E) count_faces[E]=Set{Int}() end
          push!(count_faces[E],F)
        end
      end
      @assert(all([length(it) in [1,2] for it in values(count_faces)]))
    end

    # one face should have at least 3 edges
    begin
      for (F,fe) in enumerate(ret.C[:FE])
        @assert(length(fe)>=3)

        # and should be able to find face cycles
        cell_ev=[ret.C[:EV][E] for E in fe]
        find_vcycles(cell_ev)
      end
    end

  end

  if debug_mode
    VIEWCOMPLEX(ret, explode=[1.0,1.0,1.0], show=["V","FV","Vtext"], title="arrange2d_v2 / output / final")
  end

  return ret
end

# /////////////////////////////////////////////////////////////////////
function arrange2d_v2(lar::Lar)
  return arrange2d_v2(lar.V,lar.C[:EV])
end
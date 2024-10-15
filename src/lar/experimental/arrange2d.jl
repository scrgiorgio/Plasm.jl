# /////////////////////////////////////////////////////////////////////
function arrange2d_v2(lar::Lar; debug_mode=false)

  tin = Triangulate.TriangulateIO()
  tin.pointlist = lar.V 
  tin.segmentlist = hcat(lar.C[:EV]...)  # constrained triangulation

  if debug_mode
    VIEWEDGES(tin.pointlist, tin.segmentlist, title="arrange2d_v2 / Triangulate input")
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
      tin.pointlist=ToPoints([p + rand(2) * LAR_EXPERIMENTAL_ARRANGE_PERTURBATION for p in eachcol(lar.V)])
    end
  end

  if debug_mode
    VIEWTRIANGLES(tout.pointlist, tout.trianglelist, title="arrange2d_v2 / Triangulate output")
  end

  # round vertices
  begin
    vmap=Dict()
    pointsdb=PointsDB() 
    for (P,p) in enumerate(eachcol(tout.pointlist))
      vmap[P]=add_point(pointsdb, round_vector(Vector{Float64}(p), digits=LAR_EXPERIMENTAL_ARRANGE_ROUND))
    end
  end

  # find cycles
  begin
    cycles=Cycles()
    segments  = Set(simplify_cells([Cell([a,b])   for (a,b  ) in eachcol(tout.segmentlist)]))
    triangles =     simplify_cells([Cell([a,b,c]) for (a,b,c) in eachcol(tout.trianglelist)])
    adjacent_triangles=find_adjacents_cells(triangles, 2, segments)
    groups=find_groups_of_cells(adjacent_triangles)
  
    for (A, triangle_ids) in enumerate(groups)
  
      # each group will form a face (can be holed and non-convex)
      complex_face=Cells()
      for triangle_id in triangle_ids 
        u,v,w = triangles[triangle_id]
        for (a,b) in [ [u,v], [v,w], [w,u] ]
          a,b = normalize_cell([a,b])
          if [a,b] in segments
            push!(complex_face, [a,b])
          end
        end
      end
  
      complex_face=simplify_cells(complex_face)
      for cycle in find_vcycles(complex_face)
        cycle=normalize_cycle( [[vmap[a],vmap[b]] for (a,b) in cycle]) # go to the rounded world, so something can get filtered
        if length(cycle)>=3
          push!(cycles, cycle)
        end
      end
      
    end


  end

  # build lar
  begin
    lar=Lar(get_points(pointsdb), Dict( :EV => Cells(), :FE => Cells()))
    for cycle in cycles
      fe=Cell()
      for (I,(a,b)) in enumerate(cycle)
        a,b= (I==length(cycle)) ? cycle[end] : [cycle[I][1],cycle[I+1][1]]
        push!(lar.C[:EV], [a,b]) # will simplify later
        push!(fe, length(lar.C[:EV]))
      end
      push!(lar.C[:FE], fe)
    end
  end

  lar=SIMPLIFY(lar)
  lar.C[:FV]=compute_FV(lar)
  CHECK(lar)

  if debug_mode
    VIEWCOMPLEX(lar, explode=[1.0,1.0,1.0], show=["V","FV","Vtext"], title="arrange2d_v2 / output / final")
  end

  return lar

end
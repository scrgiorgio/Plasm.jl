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
      println("WARNING Triangulate.triangulate failed, so perturnbing the points")
      tin.pointlist=ToPoints([p + rand(2) * LAR_ARRANGE2D_PERTURBATION for p in eachcol(lar.V)])
    end
  end

  if debug_mode
    VIEWTRIANGLES(tout.pointlist, tout.trianglelist, title="arrange2d_v2 / Triangulate output")
  end

  # NOTE: segment list does not contain internal edges, but only "important" edges
  lar=Lar(tout.pointlist)
  segments  = Set(simplify_cells([Cell([a,b])   for (a,b  ) in eachcol(tout.segmentlist)]))
  triangles =     simplify_cells([Cell([a,b,c]) for (a,b,c) in eachcol(tout.trianglelist)])

  # round vertices
  begin
    vmap=Dict()
    pointsdb=points_db=PointsDB() 
    for (P,p) in enumerate(eachcol(lar.V))
      vmap[P]=add_point(pointsdb, round_vector(Vector{Float64}(p), digits=LAR_ARRANGE2D_ROUND))
    end
    lar.V=get_points(pointsdb)
    segments   = Set([it for it in simplify_cells( [[vmap[a],vmap[b]]         for (a,b)   in segments ]) if length(it)==2])
    triangles  =     [it for it in simplify_cells( [[vmap[a],vmap[b],vmap[c]] for (a,b,c) in triangles]) if length(it)==3]
  end

  # build faces
  begin
    lar.C[:FE]=Cells() 
    lar.C[:EV]=Cells()
    for cycle in find_triangles_cycles(triangles, segments)
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
  COMPUTE(lar, :FV)
  CHECK(lar)

  if debug_mode
    VIEWCOMPLEX(lar, explode=[1.0,1.0,1.0], show=["V","FV","Vtext"], title="arrange2d_v2 / output / final")
  end

  return lar

end
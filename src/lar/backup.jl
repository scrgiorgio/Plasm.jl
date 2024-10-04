


# /////////////////////////////////////////////////////////////////////
"""

not working (!)

function arrange3d_v2(lar::Lar; debug_mode=false)

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
    VIEWCOMPLEX(lar, explode=[1.2,1.2,1.2], show=["V","EV","FV","Vtext", "Ftext"])
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
export arrange3d_v2#
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
  move_out = 3*LinearAlgebra.norm(b2 - b1)

  for attempt in 1:max_attempts

    inside_bbox    = [random_float(b1[I],b2[I]) for I in 1:length(b1) ]
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
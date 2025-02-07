
# /////////////////////////////////////////////////////////////////////////////
mutable struct GetOrientedTriangles  <: Viewer
  db::PointsDB
  triangles::Vector{Cell}
  function GetOrientedTriangles()
    new(PointsDB(), Vector{Cell}())
  end
end
export GetOrientedTriangles

# ///////////////////////////////////////////////////////////////////////
function render_points(viewer::GetOrientedTriangles, points::Vector{Float32}; colors=nothing, point_size=DEFAULT_POINT_SIZE, point_color=DEFAULT_POINT_COLOR)
  # not needed
end

# ///////////////////////////////////////////////////////////////////////
function render_lines(viewer::GetOrientedTriangles, points::Vector{Float32}; colors=nothing, line_width=DEFAULT_LINE_WIDTH, line_color=DEFAULT_LINE_COLOR)
  # not needed
end

# ///////////////////////////////////////////////////////////////////////
function add_triangle(viewer::GetOrientedTriangles, p1,p2,p3)
  a=add_point(viewer.db, p1)
  b=add_point(viewer.db, p2)
  c=add_point(viewer.db, p3)
  if a==b || b==c || c==a return end
  push!(viewer.triangles,[a,b,c])
end

# ///////////////////////////////////////////////////////////////////////
function render_triangles(viewer::GetOrientedTriangles, points::Vector{Float32}; normals=nothing, colors=nothing, face_color=DEFAULT_FACE_COLOR, enable_polygon_offset=false)
	if length(points)==0  return end
  # @show(points)
  for (a,b,c, d,e,f, g,h,i) in [points[i:i+8] for i in 1:9:length(points)]
    p1=PointNd([a,b,c])
    p2=PointNd([d,e,f])
    p3=PointNd([g,h,i])
    add_triangle(viewer,p1,p2,p3)
  end
end


# ///////////////////////////////////////////////////////////////////////
function run_viewer(viewer::GetOrientedTriangles; properties::Properties=Properties(), use_thread=true)

  points=get_points(viewer.db)
  triangles=viewer.triangles

  # only boundary triangles
  begin
    count=Dict{Cell,Int}()
    for triangle in triangles
      triangle=normalize_cell(triangle)
      if !haskey(count,triangle) count[triangle]=0 end
      count[triangle]+=1
    end
    triangles=[triangle for triangle in triangles if count[normalize_cell(triangle)]==1]
  end

  return find_connected_components(points, triangles)
end

# ///////////////////////////////////////////////////
function find_connected_components(points,triangles)

  # compute edge connectivity
  begin
    connections = Dict()

    for (a, b, c) in triangles
      
        # direct order
        for (A, B) in [(a, b), (b, c), (c, a)]
            if !haskey(connections, (A, B)) connections[(A, B)] = [] end
            push!(connections[(A, B)], (a, b, c))
        end

        # inverse order
        for (A, B) in [(a, c), (c, b), (b, a)]
            if !haskey(connections, (A, B)) connections[(A, B)] = [] end
            push!(connections[(A, B)], (a, c, b))
        end
    end
  end

  begin
    components = Vector{Cells}()
    visited = Set()
    for (a, b, c) in triangles
        if normalize_cell([a, b, c]) in visited continue end
        component,stack = Cells(),[(a, b, c)]
        while !isempty(stack)
            (a, b, c) = pop!(stack)
            if normalize_cell([a, b, c]) in visited continue end
            push!(component, [a, b, c])
            push!(visited, normalize_cell([a, b, c]))
            for (A, B) in [(a, b), (b, c), (c, a)]
              # need to go in the reverse order
              for (d, e, f) in get(connections, (B, A), [])
                  push!(stack, (d, e, f))
              end
            end
        end
        push!(components, component)
    end
  end

  return (points,components)
end
  
export find_connected_components

# ///////////////////////////////////////////////////
function get_oriented_triangles(hpc::Hpc)
  viewer=GetOrientedTriangles()
  render_hpc(viewer, hpc)
  return run_viewer(viewer)
end

export get_oriented_triangles

# WRONG (!)
function BREP(obj::Hpc)
  return CONS([S1,CAT∘S2])(get_oriented_triangles(obj))
end

# ///////////////////////////////////////////////////
function  BREP(lar::Lar) 
  # OLD AND WRONG
  #   return CONS([S1,CAT∘S2])(get_oriented_triangles(MKPOL(obj.V, obj.C[:CV])))

  points,triangles=lar.V,[]
  for face in lar.C[:FV]
    for I in 1:(length(face)-2)
      push!(triangles,[face[1], face[I+1], face[I+2]])
    end
  end

  return find_connected_components(points,triangles)
end

export BREP








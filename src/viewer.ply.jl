
# /////////////////////////////////////////////////////////////////////////////
mutable struct PLYViewer  <: Viewer
  filename::Any
  points::Vector{Float32}
  colors::Vector{Float32}
  function PLYViewer(filename)
    new(
      filename, 
      Vector{Float32}(), 
      Vector{Float32}()
    )
  end
end
export PLYViewer


# ///////////////////////////////////////////////////////////////////////
function render_points(viewer::PLYViewer, points::Vector{Float32}; colors=nothing, point_size=DEFAULT_POINT_SIZE, point_color=DEFAULT_POINT_COLOR)
	if length(points)==0 return end
  return # not supported
	
end

# ///////////////////////////////////////////////////////////////////////
function render_lines(viewer::PLYViewer, points::Vector{Float32}; colors=nothing, line_width=DEFAULT_LINE_WIDTH, line_color=DEFAULT_LINE_COLOR)
	if length(points)==0 return end
  return # not supported
end

# ///////////////////////////////////////////////////////////////////////
function render_triangles(viewer::PLYViewer, points::Vector{Float32}; normals=nothing, colors=nothing, face_color=DEFAULT_FACE_COLOR, enable_polygon_offset=false)
	if length(points)==0 return end
  if isnothing(colors)
    colors=Vector{Float32}() 
    for I in 1:div(length(points),3)
      append!(colors, face_color)
    end
  end
  append!(viewer.points, points)
  append!(viewer.colors, colors)
  # ignoring normals here
end

# ///////////////////////////////////////////////////////////////////////
function run_viewer(viewer::PLYViewer; properties::Properties=Properties(), use_thread=true)

  num_points    =div(length(viewer.points),3)
  num_triangles =div(num_points,3)

  file=open(viewer.filename, "w")
  write(file, "ply\n")
  write(file, "format ascii 1.0\n")
  write(file, "comment author: Plasm.JL\n")
  write(file, "element vertex $(num_points)\n")
  write(file, "property float x\n")
  write(file, "property float y\n")
  write(file, "property float z\n")
  write(file, "property uchar red\n")
  write(file, "property uchar green\n")
  write(file, "property uchar blue\n")
  write(file, "property uchar alpha\n")
  write(file, "element face $(num_triangles)\n")
  write(file, "property list uchar int vertex_indices\n")
  write(file, "end_header\n")

  C=0
  for (x,y,z) in Iterators.partition(viewer.points,3)
    r,g,b,a = [Int(it*255) for it in viewer.colors[C+1:C+4]]
    write(file, "$(x) $(y) $(z) $(r) $(g) $(b) $(a)\n") 
    C+=4
  end

  for (i,j,k) in Iterators.partition(1:num_points,3)
    write(file, "3 $(i-1) $(j-1) $(k-1)\n")
  end

  close(file)
end
export run_viewer
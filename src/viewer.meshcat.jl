using MeshCat: Visualizer
using CoordinateTransformations
using GeometryBasics: Point, TriangleFace, Mesh, meta
using Colors

# ///////////////////////////////////////////////////////////////////////
mutable struct MeshCatViewer <: Viewer
	vis::Visualizer
	objects::Vector{String}
	bbox::Box3d

	# constructor
	function MeshCatViewer(notebook=false)
		println("MeshCatViewer",notebook)
		vis=Visualizer()
		self=new(vis,Vector{String}(),invalidBox())

		if !notebook
			# render(vis) better to do in another cell
			open(vis)
			wait(vis)
		end

		return self
	end

end
export MeshCatViewer


# ///////////////////////////////////////////////////////////////////////
function render_points(viewer::MeshCatViewer, points::Vector{Float32}; colors=nothing, point_size=DEFAULT_POINT_SIZE, point_color=DEFAULT_POINT_COLOR)
	if length(points)==0 return end
	vis=viewer.vis

	for (x,y,z) in split(points,3)
		addPoint(viewer.bbox, Point3d(x,y,z))
	end

	points=[Point(x,y,z) for (x,y,z) in split(points,3)]

	# points, always use colored points otherwise all black
	if isnothing(colors)
		colors=[RGBA(point_color...) for I in 1:length(points)]
	else
		colors=[RGBA(r,g,b,a) for (r,g,b,a) in split(colors,4)]
	end

	id=string(rand(Int))
	geometry = PointCloud(points, colors)
	setobject!(vis[id], geometry,PointsMaterial(color=RGBA(1.0,1.0,1.0,1.0),size=point_size/100.0))


	push!(viewer.objects,id)

end
export render_points

# ///////////////////////////////////////////////////////////////////////
function render_lines(viewer::MeshCatViewer, points::Vector{Float32}; colors=nothing, line_width=DEFAULT_LINE_WIDTH, line_color=DEFAULT_LINE_COLOR)
	if length(points)==0 return end
	vis=viewer.vis

	for (x,y,z) in split(points,3)
		addPoint(viewer.bbox, Point3d(x,y,z))
	end

  # Due to limitations of the OpenGL Core Profile 55 with the WebGL renderer on most platforms linewidth will always be 1 regardless of the set value.
	points=Vector{Point{3, Float64}}([Point(x,y,z) for (x,y,z) in split(points,3) ])

	if isnothing(colors)
		geometry=LineSegments(points, LineBasicMaterial(color=RGBA(line_color...), vertexColors=false, linewidth =line_width/100.0))
	else
		colors=[RGBA(r,g,b,a) for (r,g,b,a) in split(colors,4)]
		geometry=LineSegments(PointCloud(points,colors), LineBasicMaterial(color=RGBA(1.0,1.0,1.0,1.0), vertexColors=true, linewidth=line_width/100.0))
	end
	
	id=string(rand(Int))
	setobject!(vis[id], geometry)
	push!(viewer.objects,id)
end
export render_lines

# ///////////////////////////////////////////////////////////////////////
function render_triangles(viewer::MeshCatViewer, points::Vector{Float32}; normals=nothing, colors=nothing, face_color=DEFAULT_FACE_COLOR, enable_polygon_offset=false)
	if length(points)==0 return end
	vis=viewer.vis

	for (x,y,z) in split(points,3)
		addPoint(viewer.bbox, Point3d(x,y,z))
	end

	id=string(rand(Int))
	if isnothing(normals)
		normals=compute_triangles_normals(points)
	end

	points =[Point(x,y,z) for (x,y,z) in split(points, 3)]
	faces  =[TriangleFace(I,I+1,I+2) for I in 1:3:length(points)]
	geometry=Mesh(points,faces)

	if isnothing(colors)
		normals=[Point(x,y,z) for (x,y,z) in split(normals,3)]
		setobject!(vis[id], meta(geometry; normals=normals), MeshLambertMaterial(color=RGBA(face_color...),vertexColors=false))
	else
		colors=[RGBA(r,g,b,a) for (r,g,b,a) in split(colors,4)]
		setobject!(vis[id], meta(geometry, normals=normals, vertexColors=colors), MeshLambertMaterial(color=RGBA(1.0,1.0,1.0,1.0),vertexColors=true))
	end

	push!(viewer.objects,id)
	
end
export render_triangles

# ///////////////////////////////////////////////////////////////////////
function run_viewer(viewer::MeshCatViewer; properties::Properties=Properties(),use_thread=false)
	
	# todo
	@assert(!use_thread)
	vis=viewer.vis

	# reduce to the case -1,+1 (do not consider the axis)
	begin
		box_size = viewer.bbox.p2 - viewer.bbox.p1
		max_size = maximum([box_size[1], box_size[2], box_size[3]])
		box_center = center(viewer.bbox)
	end

	show_axis = get(properties, "show_axis", true)
	if show_axis
		render_axis(viewer, Point3d(0.0, 0.0, 0.0), Point3d(1.05, 1.05, 1.05))
	end

	for id in viewer.objects
		settransform!(vis[id], LinearMap(UniformScaling(2.0/max_size)) ∘ Translation(-box_center...))
	end

	setvisible!(vis["/Grid"], false)
	setvisible!(vis["/Axes"], false)

	setprop!(vis["/Lights/AmbientLight/<object>"], "intensity", 0.88)

	background_color = get(properties, "background_color", DEFAULT_BACKGROUND_COLOR)
  setprop!(vis["/Background"], "top_color",    RGBA(background_color...))
  setprop!(vis["/Background"], "bottom_color", RGBA(background_color...))

	# any attempt to set the camera was unsuccessful
	# https://github.com/ferrolho/meshcat#camera-control

  setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 1.0)
  setprop!(vis["/Cameras/default/rotated/<object>"], "fov", get(properties, "fov", DEFAULT_FOV))
  setprop!(vis["/Cameras/default/rotated/<object>"], "near", 0.01)
  setprop!(vis["/Cameras/default/rotated/<object>"], "far", 10.0)

end
export run_viewer
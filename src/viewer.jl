Properties = Dict{String,Any}
export Properties


# /////////////////////////////////////////////////////////////////////
mutable struct GLVertexBuffer

	id::Int32
	vector::Vector{Float32}

  # constructor
  function GLVertexBuffer()
    ret = new(-1, Vector{Float32}())
    finalizer(releaseGpuResources, ret)
    return ret
  end

  # constructor
  function GLVertexBuffer(vector::Vector{Float32})
    ret = new(-1, copy(vector))
    finalizer(releaseGpuResources, ret)
    return ret
  end

end

export GLVertexBuffer


# /////////////////////////////////////////////////////////////////////
mutable struct GLVertexArray

	id::Int32

  # constructor
  function GLVertexArray()
    ret = new(-1)
    finalizer(releaseGpuResources, ret)
    return ret
  end

end
export GLVertexArray


# /////////////////////////////////////////////////////////////////////
mutable struct GLBatch

	primitive::UInt32
	T::Matrix4d
	vertex_array::GLVertexArray

	vertices::GLVertexBuffer
	normals::GLVertexBuffer
	colors::GLVertexBuffer

	point_size::Int
	line_width::Int
	point_color::Point4d
	line_color::Point4d
	face_color::Point4d

	enable_polygon_offset::Bool

  function GLBatch(prim::UInt32=GL_POINTS)
    ret = new(
      prim,
      Matrix4d(),
      GLVertexArray(),
      GLVertexBuffer(),
      GLVertexBuffer(),
      GLVertexBuffer(),
      1,
      1,
      Point4d(0.0, 0.0, 0.0, 1.0),
      Point4d(1.0, 1.0, 1.0, 1.0),
      Point4d(0.5, 0.5, 0.5, 1.0),
      false
    )
    finalizer(releaseGpuResources, ret)
    return ret
  end

end
export GLBatch


# ///////////////////////////////////////////////////////////////////////
function prependTransformation(self::GLBatch, T::Matrix4d)
	self.T = T * self.T
end
export prependTransformation

# ///////////////////////////////////////////////////////////////////////
function GetBoundingBox(batch::GLBatch)
	box = invalidBox()
	vertices = batch.vertices.vector
	for I in 1:3:length(vertices)
		point = Point3d(vertices[I+0], vertices[I+1], vertices[I+2])
		addPoint(box, point)
	end
	return box
end
export GetBoundingBox


# ////////////////////////////////////////////////////////////////////////
function GLCuboid(box::Box3d)
	points = getPoints(box)

	faces = [[1, 2, 3, 4], [4, 3, 7, 8], [8, 7, 6, 5], [5, 6, 2, 1], [6, 7, 3, 2], [8, 5, 1, 4]]

	vertices = Vector{Float32}()
	normals = Vector{Float32}()
	for face in faces

		p3, p2, p1, p0 = points[face[1]], points[face[2]], points[face[3]], points[face[4]] # reverse order
		n = 0.5 * (computeNormal(p0, p1, p2) + computeNormal(p0, p2, p3))

		append!(vertices, p0)
		append!(normals, n)
		append!(vertices, p1)
		append!(normals, n)
		append!(vertices, p2)
		append!(normals, n)
		append!(vertices, p0)
		append!(normals, n)
		append!(vertices, p2)
		append!(normals, n)
		append!(vertices, p3)
		append!(normals, n)
	end

	ret = GLBatch(GL_TRIANGLES)
	ret.vertices = GLVertexBuffer(vertices)
	ret.normals = GLVertexBuffer(normals)
	return ret
end
export GLCuboid 


# ////////////////////////////////////////////////////////////////////////
function GLAxis(p0::Point3d, p1::Point3d)

	vertices = Vector{Float32}()
	colors = Vector{Float32}()

	R = Point4d(1, 0, 0, 1)
	G = Point4d(0, 1, 0, 1)
	B = Point4d(0, 0, 1, 1)

	append!(vertices, p0)
	append!(vertices, Point3d(p1[1], p0[2], p0[3]))
	append!(colors, R)
	append!(colors, R)

	append!(vertices, p0)
	append!(vertices, Point3d(p0[1], p1[2], p0[3]))
	append!(colors, G)
	append!(colors, G)

	append!(vertices, p0)
	append!(vertices, Point3d(p0[1], p0[2], p1[3]))
	append!(colors, B)
	append!(colors, B)

	ret = GLBatch(GL_LINES)
	ret.vertices = GLVertexBuffer(vertices)
	ret.colors = GLVertexBuffer(colors)
	return ret
end

export GLAxis


# /////////////////////////////////////////////////////////////////////////////
mutable struct Viewer
	win::Any
	W::Int32
	H::Int32
	scalex::Float64
	scaley::Float64
	fov::Float64
	pos::Point3d
	dir::Point3d
	vup::Point3d
	zNear::Float64
	zFar::Float64
	walk_speed::Float64
	mouse_beginx::Float64
	mouse_beginy::Float64
	down_button::Int32
	batches::Any
	shaders::Dict{Any,Any}
	use_ortho::Bool
	exitNow::Bool
	show_lines::Bool
	background_color::Vector{Float64}
	title::String
	lighting_enabled::Bool

	# constructor
	function Viewer(batches)
		new(
			0,
			1024, 768,
			1.0, 1.0,
			DEFAULT_FOV,
			Point3d(0, 0, 1),
			Point3d(0, 0, -1),
			Point3d(1, 0, 0),
			0.01, 100.0, 0.1,
			0, 0, 0,
			batches,
			Dict{Any,Any}(),
			DEFAULT_USE_ORTHO,
			false,  # exitNow
			DEFAULT_SHOW_LINES,
			DEFAULT_BACKGROUND_COLOR,
			"Plasm.jl",
			DEFAULT_LIGHTING_ENABLED
		)
	end

end

# ///////////////////////////////////////////////////////////////////////
function GLView(batches::Vector{GLBatch}; properties::Properties=Properties())

	batches=[batch for batch in batches if length(batch.vertices.vector)>0]

	# calculate bounding box
	BOX::Box3d = invalidBox()
	for batch in batches
		box = GetBoundingBox(batch)
		addPoint(BOX, box.p1)
		addPoint(BOX, box.p2)
	end

	Size = BOX.p2 - BOX.p1
	MaxSize = max(Size[1], Size[2], Size[3])
	Center = center(BOX)

	show_axis = get(properties, "show_axis", true)
	if show_axis
		push!(batches, GLAxis(Point3d(0, 0, 0), Point3d(2, 2, 2)))
	end

	#@show(BOX)
	default_use_ortho=Size[3]==0.0
	default_pos=Center + 3.0*Point3d(default_use_ortho ? 0 : MaxSize, default_use_ortho ? 0 : MaxSize, MaxSize) 
	default_dir=normalized(Center - default_pos)
	default_vup=default_use_ortho ? Point3d(0, 1, 0) : Point3d(0, 0, 1)

	default_znear      = MaxSize * 0.001
	default_zfar       = MaxSize * 10.0
	default_walk_speed = MaxSize * 0.01 

	viewer = Viewer(batches)
	viewer.background_color =            get(properties, "background_color", viewer.background_color)
	viewer.title            =            get(properties, "title", viewer.title)
	viewer.use_ortho        =            get(properties, "use_ortho", default_use_ortho)
	viewer.show_lines       =            get(properties, "show_lines", viewer.show_lines)
	viewer.fov              =            get(properties, "fov", viewer.fov)
	viewer.pos              =            get(properties, "pos", default_pos)
	viewer.dir              = normalized(get(properties, "dir", default_dir))
	viewer.vup              = normalized(get(properties, "vup", default_vup))
	viewer.zNear            =            get(properties, "znear", default_znear)
	viewer.zFar             =            get(properties, "zfar ", default_zfar)
	viewer.walk_speed       =            get(properties, "walk_speed", default_walk_speed)
	viewer.lighting_enabled =            get(properties, "lighting_enabled", viewer.lighting_enabled)

	RunViewer(viewer)

end
export GLView
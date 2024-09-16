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
function prependTransformation!(T::Matrix4d, batch::GLBatch)
	# apply tranformation
	vertices = batch.vertices.vector
	for I in 1:3:length(vertices)
		x,y,z,w = vertices[I:I+2]...,1.0
		x,y,z,w=[
			T[1,1]*x + T[1,2]*y + T[1,3]*z +  T[1,4]*w,
			T[2,1]*x + T[2,2]*y + T[2,3]*z +  T[2,4]*w,
			T[3,1]*x + T[3,2]*y + T[3,3]*z +  T[3,4]*w,
			T[4,1]*x + T[4,2]*y + T[4,3]*z +  T[4,4]*w
		]
		vertices[I:I+2].=x/w,y/w,z/w
	end
end

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

using LinearAlgebra
using ModernGL
using GLFW
using StaticArrays

export Point3d, Point4d, Box3d, Matrix3d, Matrix4d, Quaternion, GLBatch,
	norm, normalized, invalidBox, addPoint, getPoints, center, flatten, translateMatrix, scaleMatrix, lookAtMatrix, perspectiveMatrix, orthoMatrix, getLookAt, GetBoundingBox,
	convertToQuaternion, convertToMatrix, prependTransformation, computeNormal, GLVertexBuffer,
	POINTS, LINES, TRIANGLES,
	GLCuboid, GLAxis, GLCells, GLView, GLExplode, COLORS, GetColorByName, get_ortho_scale, set_ortho_scale,
	Properties


# /////////////////////////////////////////////////////////////////////
__ortho_scale::Float64 = 1.0

function get_ortho_scale()
	global __ortho_scale
	return __ortho_scale
end

function set_ortho_scale(value::Float64)
	global __ortho_scale
	__ortho_scale = value
end


import Base: *
import Base.:-
import Base.:+

# //////////////////////////////////////////////////////////////////////////////
const Point3d = MVector{3,Float64}
Point3d() = Point3d(0.0, 0.0, 0.0)


function norm(p::Point3d)
	return sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3])
end

function normalized(p::Point3d)
	len = norm(p)
	return Point3d(p[1] / len, p[2] / len, p[3] / len)
end

function +(a::Point3d, b::Point3d)
	return Point3d(a[1] + b[1], a[2] + b[2], a[3] + b[3])
end

function -(a::Point3d, b::Point3d)
	return Point3d(a[1] - b[1], a[2] - b[2], a[3] - b[3])
end

function *(a::Point3d, vs::Float64)
	return Point3d(a[1] * vs, a[2] * vs, a[3] * vs)
end

# //////////////////////////////////////////////////////////////////////////////
const Point4d = MVector{4,Float64}

Point4d() = Point4d(0.0, 0.0, 0.0, 0.0)

function norm(p::Point4d)
	return sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + p[4] * p[4])
end

function normalized(p::Point4d)
	len = norm(p)
	return Point3d(p[1] / len, p[2] / len, p[3] / len, p[4] / len)
end

function dropW(p::Point4d)
	return Point3d(p[1] / p[4], p[2] / p[4], p[3] / p[4])
end

function +(a::Point4d, b::Point4d)
	return Point3d(a[1] + b[1], a[2] + b[2], a[3] + b[3], a[4] + b[4])
end

function -(a::Point4d, b::Point4d)
	return Point3d(a[1] - b[1], a[2] - b[2], a[3] - b[3], a[4] - b[4])
end

function *(a::Point4d, vs::Float64)
	return Point3d(a[1] * vs, a[2] * vs, a[3] * vs, a[4] * vs)
end


# ////////////////////////////////////////////////////////////////////////////
mutable struct Box3d
	p1::Point3d
	p2::Point3d

	# constyructor
	function Box3d()
		new(Point3d(0, 0, 0), Point3d(0, 0, 0))
	end

	# constructor
	function Box3d(p1::Point3d, p2::Point3d)
		new(p1, p2)
	end

end

Base.:(==)(a::Box3d, b::Box3d) = a.p1 == b.p1 && a.p2 == b.p2

function invalidBox()
	m, M = typemin(Float64), typemax(Float64)
	return Box3d(Point3d(M, M, M), Point3d(m, m, m))
end

function addPoint(box::Box3d, point::Point3d)
	for i in 1:3
		box.p1[i] = min(box.p1[i], point[i])
		box.p2[i] = max(box.p2[i], point[i])
	end
	return box
end

function getPoints(box::Box3d)
	return [
		Point3d(box.p1[1], box.p1[2], box.p1[3]),
		Point3d(box.p2[1], box.p1[2], box.p1[3]),
		Point3d(box.p2[1], box.p2[2], box.p1[3]),
		Point3d(box.p1[1], box.p2[2], box.p1[3]),
		Point3d(box.p1[1], box.p1[2], box.p2[3]),
		Point3d(box.p2[1], box.p1[2], box.p2[3]),
		Point3d(box.p2[1], box.p2[2], box.p2[3]),
		Point3d(box.p1[1], box.p2[2], box.p2[3])
	]
end

function center(box::Box3d)
	return Point3d(
		box.p2[1] - box.p1[1],
		box.p2[2] - box.p1[2],
		box.p2[3] - box.p1[3])
end


const Matrix3d = MMatrix{3,3,Float64}
Matrix3d() = Matrix3d([
	1.0 0.0 0.0;
	0.0 1.0 0.0;
	0.0 0.0 1.0])
Matrix3d(
	a0::Float64, a1::Float64, a2::Float64,
	a3::Float64, a4::Float64, a5::Float64,
	a6::Float64, a7::Float64, a8::Float64) = Matrix3d([a0 a1 a2; a3 a4 a5; a6 a7 a8])

function Matrix3d(v::Vector{Float32})
	@assert length(v) == 9
	return Matrix3d(v...)
end

function flatten(T::Matrix3d)
	return Vector{Float32}([
		T[1, 1], T[1, 2], T[1, 3],
		T[2, 1], T[2, 2], T[2, 3],
		T[3, 1], T[3, 2], T[3, 3]])
end




# ////////////////////////////////////////////////////////////////////
const Matrix4d = MMatrix{4,4,Float64}
Matrix4d() = Matrix4d([
	1.0 0.0 0.0 0.0;
	0.0 1.0 0.0 0.0;
	0.0 0.0 1.0 0.0;
	0.0 0.0 0.0 1.0])
Matrix4d(
	a0::Float64, a1::Float64, a2::Float64, a3::Float64,
	a4::Float64, a5::Float64, a6::Float64, a7::Float64,
	a8::Float64, a9::Float64, a10::Float64, a11::Float64,
	a12::Float64, a13::Float64, a14::Float64, a15::Float64) = Matrix4d([a0 a1 a2 a3; a4 a5 a6 a7; a8 a9 a10 a11; a12 a13 a14 a15])


function Matrix4d(v::Vector{Float64})
	@assert length(v) == 16
	return Matrix4d(v...)
end

function dropW(T::Matrix4d)
	return Matrix3d(
		T[1, 1], T[1, 2], T[1, 3],
		T[2, 1], T[2, 2], T[2, 3],
		T[3, 1], T[3, 2], T[3, 3])
end

function flatten(T::Matrix4d)
	return Vector{Float32}([
		T[1, 1], T[1, 2], T[1, 3], T[1, 4],
		T[2, 1], T[2, 2], T[2, 3], T[2, 4],
		T[3, 1], T[3, 2], T[3, 3], T[3, 4],
		T[4, 1], T[4, 2], T[4, 3], T[4, 4]])
end

function translateMatrix(vt::Point3d)
	return Matrix4d(
		1.0, 0.0, 0.0, vt[1],
		0.0, 1.0, 0.0, vt[2],
		0.0, 0.0, 1.0, vt[3],
		0.0, 0.0, 0.0, 1.0)
end

function scaleMatrix(vs::Point3d)
	return Matrix4d(
		vs[1], 0.0, 0.0, 0.0,
		0.0, vs[2], 0.0, 0.0,
		0.0, 0.0, vs[3], 0.0,
		0.0, 0.0, 0.0, 1.0)
end

function lookAtMatrix(eye::Point3d, center::Point3d, up::Point3d)
	forward = normalized(center - eye)
	side = normalized(cross(forward, up))
	up = cross(side, forward)
	m = Matrix4d(
		side[1], up[1], -forward[1], 0.0,
		side[2], up[2], -forward[2], 0.0,
		side[3], up[3], -forward[3], 0.0,
		0.0, 0.0, 0.0, 1.0
	)
	return transpose(m) * translateMatrix(-1 * eye)
end


function perspectiveMatrix(fovy::Float64, aspect::Float64, zNear::Float64, zFar::Float64)
	radians = deg2rad(fovy / 2.0)
	cotangent = cos(radians) / sin(radians)
	m = Matrix4d()
	m[1, 1] = cotangent / aspect
	m[2, 2] = cotangent
	m[3, 3] = -(zFar + zNear) / (zFar - zNear)
	m[3, 4] = -1.0
	m[4, 3] = -2.0 * zNear * zFar / (zFar - zNear)
	m[4, 4] = 0.0
	return transpose(m)
end


function orthoMatrix(left::Float64, right::Float64, bottom::Float64, top::Float64, nearZ::Float64, farZ::Float64)
	m = Matrix4d()
	m[1, 1] = 2 / (right - left)
	m[1, 2] = 0
	m[1, 3] = 0
	m[1, 4] = -(right + left) / (right - left)
	m[2, 1] = 0
	m[2, 2] = 2 / (top - bottom)
	m[2, 3] = 0
	m[2, 4] = -(top + bottom) / (top - bottom)
	m[3, 1] = 0
	m[3, 2] = 0
	m[3, 3] = -2 / (farZ - nearZ)
	m[3, 4] = -(farZ + nearZ) / (farZ - nearZ)
	m[4, 1] = 0
	m[4, 2] = 0
	m[4, 3] = 0
	m[4, 4] = 1
	return m
end

function getLookAt(T::Matrix4d)
	T = inv(T)
	pos = (Point3d(T[1, 4], T[2, 4], T[3, 4]))
	dir = normalized(Point3d(-T[1, 3], -T[2, 3], -T[3, 3]))
	vup = normalized(Point3d(T[1, 2], T[2, 2], T[3, 2]))
	pos, dir, vup
end

# //////////////////////////////////////////////////////////////////////////////
# see https://github.com/JuliaGeometry/Quaternions.jl/blob/master/src/Quaternion.jl
struct Quaternion

	w::Float64
	x::Float64
	y::Float64
	z::Float64

	# constructor
	function Quaternion()
		new(1, 0, 0, 0)
	end

	# constructor
	function Quaternion(w::Float64, x::Float64, y::Float64, z::Float64)
		norm = sqrt(w * w + x * x + y * y + z * z)
		w /= norm
		x /= norm
		y /= norm
		z /= norm
		new(w, x, y, z)
	end

	# constructor
	function Quaternion(axis::Point3d, angle::Float64)
		axis = normalized(axis)
		c = cos(0.5 * angle)
		s = sin(0.5 * angle)
		new(c, s * axis[1], s * axis[2], s * axis[3])
	end

end

(*)(q::Quaternion, w::Quaternion) = Quaternion(
	q.w * w.w - q.x * w.x - q.y * w.y - q.z * w.z,
	q.w * w.x + q.x * w.w + q.y * w.z - q.z * w.y,
	q.w * w.y - q.x * w.z + q.y * w.w + q.z * w.x,
	q.w * w.z + q.x * w.y - q.y * w.x + q.z * w.w)


function convertToQuaternion(T::Matrix4d)
	# See https://arxiv.org/pdf/math/0701759.pdf
	a2 = 1.0 + T[1, 1] + T[2, 2] + T[3, 3]
	b2 = 1.0 + T[1, 1] - T[2, 2] - T[3, 3]
	c2 = 1.0 - T[1, 1] + T[2, 2] - T[3, 3]
	d2 = 1.0 - T[1, 1] - T[2, 2] + T[3, 3]

	if a2 >= max(b2, c2, d2)
		a = sqrt(a2) / 2
		return Quaternion(a, (T[3, 2] - T[2, 3]) / 4a, (T[1, 3] - T[3, 1]) / 4a, (T[2, 1] - T[1, 2]) / 4a)
	elseif b2 >= max(c2, d2)
		b = sqrt(b2) / 2
		return Quaternion((T[3, 2] - T[2, 3]) / 4b, b, (T[2, 1] + T[1, 2]) / 4b, (T[1, 3] + T[3, 1]) / 4b)
	elseif c2 >= d2
		c = sqrt(c2) / 2
		return Quaternion((T[1, 3] - T[3, 1]) / 4c, (T[2, 1] + T[1, 2]) / 4c, c, (T[3, 2] + T[2, 3]) / 4c)
	else
		d = sqrt(d2) / 2
		return Quaternion((T[2, 1] - T[1, 2]) / 4d, (T[1, 3] + T[3, 1]) / 4d, (T[3, 2] + T[2, 3]) / 4d, d)
	end
end

function convertToMatrix(q::Quaternion)
	sx, sy, sz = 2q.w * q.x, 2q.w * q.y, 2q.w * q.z
	xx, xy, xz = 2q.x^2, 2q.x * q.y, 2q.x * q.z
	yy, yz, zz = 2q.y^2, 2q.y * q.z, 2q.z^2
	return Matrix4d(
		1 - (yy + zz), xy - sz, xz + sy, 0.0,
		xy + sz, 1 - (xx + zz), yz - sx, 0.0,
		xz - sy, yz + sx, 1 - (xx + yy), 0.0,
		0.0, 0.0, 0.0, 1.0)
end

# /////////////////////////////////////////////////////////////////////
function glGenBuffer()
	id = GLuint[0]
	glGenBuffers(1, id)
	glCheckError("generating a buffer, array, or texture")
	id[]
end

# /////////////////////////////////////////////////////////////////////
function glGenVertexArray()
	id = GLuint[0]
	glGenVertexArrays(1, id)
	glCheckError("generating a buffer, array, or texture")
	id[]
end

# /////////////////////////////////////////////////////////////////////
function glCheckError(actionName="")
	message = glErrorMessage()
	if length(message) > 0
		if length(actionName) > 0
			error("Error ", actionName, ": ", message)
		else
			error("Error: ", message)
		end
	end
end

# /////////////////////////////////////////////////////////////////////
function glErrorMessage()
	err = glGetError()
	if err == GL_NO_ERROR
		return ""
	end
	if err == GL_INVALID_ENUM
		return "GL_INVALID_ENUM"
	end
	if err == GL_INVALID_VALUE
		return "GL_INVALID_VALUE"
	end
	if err == GL_INVALID_OPERATION
		return "GL_INVALID_OPERATION"
	end
	if err == GL_INVALID_FRAMEBUFFER_OPERATION
		return "GL_INVALID_FRAMEBUFFER_OPERATION"
	end
	if err == GL_OUT_OF_MEMORY
		return "GL_OUT_OF_MEMORY"
	end
	return "Unknown OpenGL error with error code"
end


# /////////////////////////////////////////////////////////////////////
__release_gpu_resources__ = []

function glDeleteLater(fun::Function)
	global __release_gpu_resources__
	append!(__release_gpu_resources__, [fun])
end

# /////////////////////////////////////////////////////////////////////
function glDeleteNow()
	global __release_gpu_resources__
	for fun in __release_gpu_resources__
		fun()
	end
end

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

# /////////////////////////////////////////////////////////////////////
function releaseGpuResources(buffer::GLVertexBuffer)
	global __release_gpu_resources__
	if buffer.id >= 0
		id = buffer.id
		buffer.id = -1
		glDeleteLater(function ()
			glDeleteBuffers(1, [id])
		end)
	end
end

# /////////////////////////////////////////////////////////////////////
function enableAttribute(location::Int32, buffer::GLVertexBuffer, num_components::Int64)
	if length(buffer.vector) == 00 || location < 0
		return
	end
	if buffer.id < 0
		buffer.id = glGenBuffer()
	end
	glBindBuffer(GL_ARRAY_BUFFER, buffer.id)
	glBufferData(GL_ARRAY_BUFFER, sizeof(buffer.vector), buffer.vector, GL_STATIC_DRAW)
	glVertexAttribPointer(location, num_components, GL_FLOAT, false, 0, C_NULL)
	glEnableVertexAttribArray(location)
	glBindBuffer(GL_ARRAY_BUFFER, 0)
end

# /////////////////////////////////////////////////////////////////////
function disableAttribute(location::Int32, buffer::GLVertexBuffer)
	if length(buffer.vector) == 00 || location < 0
		return
	end
	glDisableVertexAttribArray(location)
end


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


# /////////////////////////////////////////////////////////////////////
function releaseGpuResources(array::GLVertexArray)
	global __release_gpu_resources__
	if array.id >= 0
		id = array.id
		array.id = -1
		glDeleteLater(function ()
			glDeleteVertexArrays(1, [id])
		end)
	end
end

# /////////////////////////////////////////////////////////////////////
function enableVertexArray(array::GLVertexArray)

	# not needed or osx
	if Sys.isapple()
		return
	end

	if array.id < 0
		array.id = glGenVertexArray()
	end
	glBindVertexArray(array.id)
end

# /////////////////////////////////////////////////////////////////////
function disableVertexArray(array::GLVertexArray)

	# not needed or osx
	if Sys.isapple()
		return
	end

	glBindVertexArray(0)
end

POINTS = GL_POINTS
LINES = GL_LINES
TRIANGLES = GL_TRIANGLES

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

	# constructor
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

function prependTransformation(self::GLBatch, T::Matrix4d)
	self.T = T * self.T
end


function releaseGpuResources(batch::GLBatch)
	releaseGpuResources(batch.vertex_array)
	releaseGpuResources(batch.vertices)
	releaseGpuResources(batch.normals)
	releaseGpuResources(batch.colors)
end



# ///////////////////////////////////////////////////////////////////////
function computeNormal(p0::Point3d, p1::Point3d, p2::Point3d)
	return normalized(cross(p1 - p0, p2 - p0))
end


# ///////////////////////////////////////////////////////////////////////
function GetBoundingBox(batch::GLBatch)
	box = invalidBox()
	vertices = batch.vertices.vector
	for I in 1:3:length(vertices)
		point = Point3d(vertices[I+0], vertices[I+1], vertices[I+1])
		addPoint(box, point)
	end
	return box
end


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

# ////////////////////////////////////////////////////////////////////////
function GLAxis(p0::Point3d, p1::Point3d)

	vertices = Vector{Float32}()
	colors = Vector{Float32}()

	R = Point4d(1, 0, 0, 1)
	append!(vertices, p0)
	append!(vertices, Point3d(p1[1], p0[2], p0[3]))
	append!(colors, R)
	append!(colors, R)
	G = Point4d(0, 1, 0, 1)
	append!(vertices, p0)
	append!(vertices, Point3d(p0[1], p1[2], p0[3]))
	append!(colors, G)
	append!(colors, G)
	B = Point4d(0, 0, 1, 1)
	append!(vertices, p0)
	append!(vertices, Point3d(p0[1], p0[2], p1[3]))
	append!(colors, B)
	append!(colors, B)

	ret = GLBatch(GL_LINES)
	ret.vertices = GLVertexBuffer(vertices)
	ret.colors = GLVertexBuffer(colors)
	return ret
end



# /////////////////////////////////////////////////////////////////////
mutable struct GLShader

	vertex_source
	frag_source

	program_id::Int32
	vertex_shader_id::Int32
	frag_shader_id::Int32

	# constructor
	function GLShader(vertex, fragment)
		ret = new(vertex, fragment, -1, -1, -1)
		finalizer(releaseGpuResources, ret)
		return ret
	end

end


# /////////////////////////////////////////////////////////////////////
function releaseGpuResources(shader::GLShader)

	global __release_gpu_resources__

	if shader.vertex_shader_id >= 0
		id = shader.vertex_shader_id
		shader.vertex_shader_id = -1
		glDeleteLater(function ()
			glDeleteShader(id)
		end)
	end

	if shader.frag_shader_id >= 0
		id = shader.frag_shader_id
		shader.frag_shader_id = -1
		glDeleteLater(function ()
			glDeleteShader(id)
		end)

	end

	if shader.program_id >= 0
		id = shader.program_id
		shader.program_id = -1
		glDeleteLater(function ()
			glDeleteProgram(id)
		end)
	end

end


# /////////////////////////////////////////////////////////////////////
function createShader(type, source)
	shader_id = glCreateShader(type)::GLuint
	glCheckError()
	glShaderSource(shader_id, 1, convert(Ptr{UInt8}, pointer([convert(Ptr{GLchar}, pointer(source))])), C_NULL)
	glCompileShader(shader_id)
	status = GLint[0]
	glGetShaderiv(shader_id, GL_COMPILE_STATUS, status)
	if status[1] == GL_FALSE
		maxlength = 8192
		buffer = zeros(GLchar, maxlength)
		sizei = GLsizei[0]
		glGetShaderInfoLog(shader_id, maxlength, sizei, buffer)
		len = sizei[]
		error_msg = unsafe_string(pointer(buffer), len)
		error("shader compilation failed\n", error_msg, "\nsource\n", source)
	end
	return shader_id
end



# /////////////////////////////////////////////////////////////////////
function enableProgram(shader)

	if (shader.program_id < 0)

		shader.program_id = glCreateProgram()
		glCheckError()

		shader.vertex_shader_id = createShader(GL_VERTEX_SHADER, shader.vertex_source)
		glAttachShader(shader.program_id, shader.vertex_shader_id)
		glCheckError()

		shader.frag_shader_id = createShader(GL_FRAGMENT_SHADER, shader.frag_source)
		glAttachShader(shader.program_id, shader.frag_shader_id)
		glCheckError()

		glLinkProgram(shader.program_id)
		glCheckError()
		status = GLint[0]
		glGetProgramiv(shader.program_id, GL_LINK_STATUS, status)
		if status[1] == GL_FALSE
			maxlength = 8192
			buffer = zeros(GLchar, maxlength)
			sizei = GLsizei[0]
			glGetProgramInfoLog(shader.program_id, maxlength, sizei, buffer)
			len = sizei[]
			error_msg = unsafe_string(pointer(buffer), len)
			error("Error linking program\n", error_msg)
		end
		glCheckError()
	end

	glUseProgram(shader.program_id)
end

# /////////////////////////////////////////////////////////////////////
function disableProgram(shader)
	glUseProgram(0)
end

# /////////////////////////////////////////////////////////////////////
include("./viewer.shader.jl")


# /////////////////////////////////////////////////////////////////////
function GLPhongShader(lighting_enabled, color_attribute_enabled)

	v = vert_source
	f = frag_source

	v = replace(v, "arg(LIGHTING_ENABLED)" => lighting_enabled ? "1" : "0")
	f = replace(f, "arg(LIGHTING_ENABLED)" => lighting_enabled ? "1" : "0")

	v = replace(v, "arg(COLOR_ATTRIBUTE_ENABLED)" => color_attribute_enabled ? "1" : "0")
	f = replace(f, "arg(COLOR_ATTRIBUTE_ENABLED)" => color_attribute_enabled ? "1" : "0")

	# this is needed for #version 330
	# v=replace(v, "attribute"=>"in")
	# f=replace(f, "attribute"=>"in")

	# v=replace(v, "varying"=>"out")
	# f=replace(f, "varying"=>"out")

	# v=string("#version 120\n",v)
	# f=string("#version 120\n",f)

	return GLShader(v, f)
end


# /////////////////////////////////////////////////////////////
mutable struct FrustumMap

	viewport::Matrix4d
	projection::Matrix4d
	modelview::Matrix4d

	inv_viewport::Matrix4d
	inv_projection::Matrix4d
	inv_modelview::Matrix4d

	# constructor
	function FrustumMap(viewport, projection::Matrix4d, modelview::Matrix4d)
		x = viewport[1]
		y = viewport[2]
		w = viewport[3]
		h = viewport[4]
		viewport_T = Matrix4d(
			w / 2.0, 0.0, 0.0, x + w / 2.0,
			0.0, h / 2.0, 0.0, y + h / 2.0,
			0.0, 0.0, 1 / 2.0, 1 / 2.0,
			0.0, 0.0, 0.0, 1.0)
		new(viewport_T, projection, modelview, inv(viewport_T), inv(projection), inv(modelview))
	end

end

function projectPoint(map::FrustumMap, p3::Point3d)
	p4 = (map.viewport * (map.projection * (map.modelview * Point4d(p3[1], p3[2], p3[3], 1.0))))
	return Point3d(p4[1] / p4[4], p4[2] / p4[4], p4[3] / p4[4])
end

function unprojectPoint(map::FrustumMap, x::Float64, y::Float64, z::Float64)
	p4 = (map.inv_modelview * (map.inv_projection * (map.inv_viewport * Point4d(x, y, z, 1.0))))
	return Point3d(p4[1] / p4[4], p4[2] / p4[4], p4[3] / p4[4])
end


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
function releaseGpuResources(viewer::Viewer)

	for batch in viewer.batches
		releaseGpuResources(batch)
	end

	for (key, shader) in viewer.shaders
		releaseGpuResources(shader)
	end
end

# ///////////////////////////////////////////////////////////////////////
function runViewer(viewer::Viewer)

	redisplay(viewer)

	ret_code = GLFW.Init()
	# println("GLFW init returned ",ret_code)

	# seems not to be needed for julia 1.x
	#GLFW.WindowHint(GLFW.CONTEXT_VERSION_MAJOR, 3)
	#GLFW.WindowHint(GLFW.CONTEXT_VERSION_MINOR, 2)
	#GLFW.WindowHint(GLFW.OPENGL_FORWARD_COMPAT, GL_TRUE)
	#GLFW.WindowHint(GLFW.OPENGL_PROFILE, GLFW.OPENGL_CORE_PROFILE)

	win = GLFW.CreateWindow(viewer.W, viewer.H, viewer.title)
	viewer.win = win
	viewer.exitNow = false
	GLFW.MakeContextCurrent(win)

	#println("Running viewer",viewer.title)
	#println("GL_SHADING_LANGUAGE_VERSION ",unsafe_string(glGetString(GL_SHADING_LANGUAGE_VERSION)))
	#println("GL_VERSION                  ",unsafe_string(glGetString(GL_VERSION)))
	#println("GL_VENDOR                   ",unsafe_string(glGetString(GL_VENDOR)))
	#println("GL_RENDERER                 ",unsafe_string(glGetString(GL_RENDERER)))

	# problem of retina
	window_size = GLFW.GetWindowSize(viewer.win)
	framebuffer_size = GLFW.GetFramebufferSize(viewer.win)
	viewer.scalex = framebuffer_size[1] / Float64(window_size[1])
	viewer.scaley = framebuffer_size[2] / Float64(window_size[2])

	GLFW.SetWindowSizeCallback(win, function (win::GLFW.Window, width::Integer, height::Integer)
		handleResizeEvent(viewer)
	end)
	GLFW.SetKeyCallback(win, function (win::GLFW.Window, key, scancode, action, mods)
		handleKeyPressEvent(viewer, key, scancode, action, mods)
	end)
	GLFW.SetCursorPosCallback(win, function (win::GLFW.Window, x, y)
		handleMouseMoveEvent(viewer, x, y)
	end)
	GLFW.SetMouseButtonCallback(win, function (win::GLFW.Window, button, action, mods)
		handleMouseButtonEvent(viewer, button, action, mods)
	end)
	GLFW.SetScrollCallback(win, function (win::GLFW.Window, dx, dy)
		handleMouseWheelEvent(viewer, dy)
	end)

	handleResizeEvent(viewer::Viewer)
	while !viewer.exitNow && !GLFW.WindowShouldClose(win)
		glRender(viewer)
		GLFW.SwapBuffers(win)
		GLFW.PollEvents()
	end

	releaseGpuResources(viewer)
	glDeleteNow()
	GLFW.DestroyWindow(win)
	GLFW.Terminate()
end

Properties = Dict{String,Any}

# ///////////////////////////////////////////////////////////////////////
function GLView(batches::Vector{GLBatch}; properties::Properties=Properties(), title::String="")

	if title!=""
		properties["title"]=title
	end

	batches=[batch for batch in batches if length(batch.vertices.vector)>0]

	show_axis = get(properties, "show_axis", true)
	if show_axis
		push!(batches, GLAxis(Point3d(0, 0, 0), Point3d(2, 2, 2)))
	end

	global viewer
	viewer = Viewer(batches)

	# calculate bounding box
	BOX = invalidBox()
	for batch in viewer.batches
		box = GetBoundingBox(batch)
		addPoint(BOX, box.p1)
		addPoint(BOX, box.p2)
	end

	Size = BOX.p2 - BOX.p1
	MaxSize = max(Size[1], Size[2], Size[3])
	Center = center(BOX)

	# backward compatible
	function _get(properties, key, default)
		return get(properties, key, get(properties, replace(key, "_" => "-"), default))
	end

	viewer.background_color = _get(properties, "background_color", viewer.background_color)

	viewer.title = _get(properties, "title", viewer.title)
	viewer.use_ortho = _get(properties, "use_ortho", viewer.use_ortho)
	viewer.show_lines = _get(properties, "show_lines", viewer.show_lines)
	viewer.fov = _get(properties, "fov", viewer.fov)
	viewer.pos = _get(properties, "pos", Center + Point3d(MaxSize, MaxSize, MaxSize) * 3.0)
	viewer.dir = normalized(_get(properties, "dir", Center - viewer.pos))
	viewer.vup = normalized(_get(properties, "vup", Point3d(0, 0, 1)))
	viewer.zNear = _get(properties, "znear", MaxSize / 50.0)
	viewer.zFar = _get(properties, "zfar ", MaxSize * 10.0)
	viewer.walk_speed = _get(properties, "walk_speed", MaxSize / 100.0)
	viewer.lighting_enabled = _get(properties, "lighting_enabled", viewer.lighting_enabled)

	runViewer(viewer)

end


# ///////////////////////////////////////////////////////////////////////
function getModelview(viewer::Viewer)
	return lookAtMatrix(viewer.pos, viewer.pos + viewer.dir, viewer.vup)
end

function getProjection(viewer::Viewer)
	ratio = viewer.W / float(viewer.H)
	if viewer.use_ortho
		# euristic that seem to work well
		Z = viewer.zNear + 0.5 * (viewer.zFar - viewer.zNear)
		right = Z * tan(deg2rad(viewer.fov / 2.0))
		left = -right
		vs = get_ortho_scale()
		# println("ortho_scale ",vs)
		top = -0.5 * (right - left) / ratio
		bottom = +0.5 * (right - left) / ratio
		return orthoMatrix(vs * left, vs * right, vs * top, vs * bottom, viewer.zNear, viewer.zFar)
	else
		return perspectiveMatrix(viewer.fov, ratio, viewer.zNear, viewer.zFar)
	end

end

# ///////////////////////////////////////////////////////////////////////
function projectPoint(viewer::Viewer, pos::Point3d)
	viewport = [0, 0, viewer.W, viewer.H]
	projection = getProjection(viewer)
	modelview = getModelview(viewer)
	map = FrustumMap(viewport, projection, modelview)
	return projectPoint(map, pos)
end

# ///////////////////////////////////////////////////////////////////////
function unprojectPoint(viewer::Viewer, x::Float64, y::Float64)
	viewport = [0, 0, viewer.W, viewer.H]
	projection = getProjection(viewer)
	modelview = getModelview(viewer)
	map = FrustumMap(viewport, projection, modelview)
	P1 = unprojectPoint(map, x, viewer.H - y, -1.0)
	P2 = unprojectPoint(map, x, viewer.H - y, +1.0)
	return normalized(P2 - P1)
end

# ///////////////////////////////////////////////////////////////////////
function getShader(viewer::Viewer, lighting_enabled, color_attribute_enabled)

	lighting_enabled = lighting_enabled && viewer.lighting_enabled

	key = (lighting_enabled, color_attribute_enabled)

	if haskey(viewer.shaders, key)
		return viewer.shaders[key]
	end

	ret = GLPhongShader(lighting_enabled, color_attribute_enabled)
	viewer.shaders[key] = ret
	return ret
end

# ///////////////////////////////////////////////////////////////////////
function glRender(viewer::Viewer)
	glEnable(GL_DEPTH_TEST)
	glDepthFunc(GL_LEQUAL)
	glDisable(GL_CULL_FACE)
	glClearDepth(1.0)
	glClearColor(viewer.background_color[1], viewer.background_color[2], viewer.background_color[3], 0.00)
	glPolygonOffset(-1.0, -1.0)

	glViewport(0, 0, viewer.W, viewer.H)
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

	PROJECTION = getProjection(viewer)
	MODELVIEW = getModelview(viewer)

	lightpos = MODELVIEW * Point4d(viewer.pos[1], viewer.pos[2], viewer.pos[3], 1.0)

	for batch in viewer.batches

		# show points
		if batch.primitive == GL_POINTS

			if batch.point_color[4] > 0.0
				glPointSize(batch.point_size)
				glRenderBatch(viewer, batch, batch.point_color, PROJECTION, MODELVIEW, lightpos)
				glPointSize(1)
			end

			# show lines
		elseif batch.primitive == GL_LINES || batch.primitive == GL_LINE_STRIP || batch.primitive == GL_LINE_LOOP

			if batch.line_color[4] > 0.0
				glLineWidth(batch.line_width)
				glRenderBatch(viewer, batch, batch.line_color, PROJECTION, MODELVIEW, lightpos)
				glLineWidth(1)
			end

			# show triangles
		else

			if viewer.show_lines && batch.line_width > 0

				# https://www.glprogramming.com/red/chapter06.html

				if batch.face_color[4] > 0.0
					glEnable(GL_POLYGON_OFFSET_FILL)
					glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
					glPolygonOffset(1.0, 1.0)
					glRenderBatch(viewer, batch, batch.face_color, PROJECTION, MODELVIEW, lightpos)
					glDisable(GL_POLYGON_OFFSET_FILL)
				end

				if batch.line_color[4] > 0.0
					glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
					glLineWidth(batch.line_width)
					glRenderBatch(viewer, batch, batch.line_color, PROJECTION, MODELVIEW, lightpos)
					glLineWidth(1)
					glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
				end

			else
				if batch.face_color[4] > 0.0
					
					glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)

					if !isnothing(batch.enable_polygon_offset)
						glEnable(GL_POLYGON_OFFSET_FILL)
						glPolygonOffset(1.0, 1.0)
						glRenderBatch(viewer, batch, batch.face_color, PROJECTION, MODELVIEW, lightpos)
					else
						glRenderBatch(viewer, batch, batch.face_color, PROJECTION, MODELVIEW, lightpos)
					end

				end
			end

		end
	end

	glCheckError()

end


function glRenderBatch(viewer::Viewer, batch::GLBatch, color, PROJECTION, MODELVIEW, lightpos)

	lighting_enabled = length(batch.normals.vector) > 0 && viewer.lighting_enabled
	color_attribute_enabled = length(batch.colors.vector) > 0

	shader = getShader(viewer, lighting_enabled, color_attribute_enabled)

	enableProgram(shader)

	projection = PROJECTION
	modelview = MODELVIEW * batch.T
	normal_matrix = dropW(transpose(inv(modelview)))

	glUniformMatrix4fv(glGetUniformLocation(shader.program_id, "u_modelview_matrix"), 1, GL_TRUE, flatten(modelview))
	glUniformMatrix4fv(glGetUniformLocation(shader.program_id, "u_projection_matrix"), 1, GL_TRUE, flatten(projection))
	glUniformMatrix3fv(glGetUniformLocation(shader.program_id, "u_normal_matrix"), 1, GL_TRUE, flatten(normal_matrix))

	u_light_position = glGetUniformLocation(shader.program_id, "u_light_position")
	if u_light_position >= 0
		glUniform3f(u_light_position, lightpos[1] / lightpos[4], lightpos[2] / lightpos[4], lightpos[3] / lightpos[4])
	end

	u_color = glGetUniformLocation(shader.program_id, "u_color")
	# println(u_color," ", color, " lighting_enabled ",lighting_enabled, " color_attribute_enabled ",color_attribute_enabled)
	if u_color >= 0 && !isnothing(color)
		glUniform4f(u_color, color[1], color[2], color[3], color[4])
	end

	enableVertexArray(batch.vertex_array)

	a_position = glGetAttribLocation(shader.program_id, "a_position")
	a_normal = glGetAttribLocation(shader.program_id, "a_normal")
	a_color = glGetAttribLocation(shader.program_id, "a_color")

	enableAttribute(a_position, batch.vertices, 3)
	enableAttribute(a_normal, batch.normals, 3)
	enableAttribute(a_color, batch.colors, 4)

	@assert length(batch.vertices.vector) % 3 == 0
	glDrawArrays(batch.primitive, 0, Int64(length(batch.vertices.vector) / 3))

	disableAttribute(a_position, batch.vertices)
	disableAttribute(a_normal, batch.normals)
	disableAttribute(a_color, batch.colors)
	disableVertexArray(batch.vertex_array)
	disableProgram(shader)

	glDepthMask(true)
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
	glDisable(GL_POLYGON_OFFSET_LINE)

end

# ///////////////////////////////////////////////////////////////////////
function redisplay(viewer::Viewer)
	# nothing to do
end

# ///////////////////////////////////////////////////////////////////////
function handleResizeEvent(viewer)
	size = GLFW.GetWindowSize(viewer.win)
	viewer.W = size[1] * viewer.scalex
	viewer.H = size[2] * viewer.scaley
	redisplay(viewer)
end

# ///////////////////////////////////////////////////////////////////////
function handleMouseButtonEvent(viewer::Viewer, button, action, mods)

	button = Dict(GLFW.MOUSE_BUTTON_1 => 1, GLFW.MOUSE_BUTTON_2 => 3, GLFW.MOUSE_BUTTON_3 => 2)[button]

	if action == GLFW.PRESS && viewer.down_button == 0
		viewer.down_button = button
		redisplay(viewer)
		return
	end

	if action == GLFW.RELEASE && button == viewer.down_button
		viewer.down_button = 0
		redisplay(viewer)
		return
	end
end

# ///////////////////////////////////////////////////////////////////////
function handleMouseMoveEvent(viewer::Viewer, x::Float64, y::Float64)

	x = x * viewer.scalex
	y = y * viewer.scaley

	button = viewer.down_button

	if (button == 0)
		viewer.mouse_beginx = x
		viewer.mouse_beginy = y
		return
	end

	deltax = float(x - viewer.mouse_beginx)
	deltay = float(viewer.mouse_beginy - y)
	W = viewer.W
	H = viewer.H

	modelview = getModelview(viewer)

	if button == 1
		screen_center = Point3d(W / 2.0, H / 2.0, 0.0)
		a = (Point3d((float)(viewer.mouse_beginx - screen_center[1]), (float)(H - viewer.mouse_beginy - screen_center[2]), 0)) * (1.0 / min(W, H))
		b = (Point3d((float)(x - screen_center[1]), (float)(H - y - screen_center[2]), 0)) * (1.0 / min(W, H))
		a[3] = 2.0^(-0.5 * norm(a))
		b[3] = 2.0^(-0.5 * norm(b))
		a = normalized(a)
		b = normalized(b)
		axis = normalized(cross(a, b))
		angle = acos(dot(a, b))

		#vt=Point3d(modelview[1,4],modelview[2,4],modelview[3,4])
		#modelview=translateMatrix(vt) * convertToMatrix(convertToQuaternion(modelview))

		q = Quaternion(axis, angle) * convertToQuaternion(modelview)
		vt = Point3d(modelview[1, 4], modelview[2, 4], modelview[3, 4])
		modelview = translateMatrix(vt) * convertToMatrix(q)

	elseif button == 3
		vt = Point3d(deltax * viewer.walk_speed, deltay * viewer.walk_speed, 0.0)
		modelview = translateMatrix(vt) * modelview
	end

	viewer.pos, viewer.dir, viewer.vup = getLookAt(modelview)

	viewer.mouse_beginx = x
	viewer.mouse_beginy = y
	redisplay(viewer)
end

# ///////////////////////////////////////////////////////////////////////
function handleMouseWheelEvent(viewer::Viewer, delta)
	viewer.pos = viewer.pos + viewer.dir * ((delta >= 0 ? 10.0 : -10.0) * viewer.walk_speed)
	redisplay(viewer)
end

# ///////////////////////////////////////////////////////////////////////
function handleKeyPressEvent(viewer::Viewer, key, scancode, action, mods)

	if action != GLFW.PRESS && action != GLFW.REPEAT
		return
	end

	if key == GLFW.KEY_ESCAPE
		viewer.exitNow = true
		return
	end

	if (key == GLFW.KEY_KP_ADD)
		viewer.walk_speed *= 0.95
		return
	end

	if (key == GLFW.KEY_KP_SUBTRACT)
		viewer.walk_speed *= (1.0 / 0.95)
		return
	end

	if (key == GLFW.KEY_W)
		dir = unprojectPoint(viewer, 0.5 * viewer.W, 0.5 * viewer.H)
		#println("dir",dir,"walk_speed",viewer.walk_speed)
		viewer.pos = viewer.pos + dir * viewer.walk_speed
		redisplay(viewer)
		return
	end

	if (key == GLFW.KEY_S)
		dir = unprojectPoint(viewer, 0.5 * viewer.W, 0.5 * viewer.H)
		viewer.pos = viewer.pos - dir * viewer.walk_speed
		redisplay(viewer)
		return
	end

	if (key == GLFW.KEY_O)
		viewer.use_ortho = !viewer.use_ortho
		#println("use_ortho ",viewer.use_ortho)
		redisplay(viewer)
		return
	end


	if (key == GLFW.KEY_UP)
		viewer.pos = viewer.pos + viewer.vup * viewer.walk_speed
		redisplay(viewer)
		return
	end

	if (key == GLFW.KEY_DOWN)
		viewer.pos = viewer.pos - viewer.vup * viewer.walk_speed
		redisplay(viewer)
		return
	end

	if (key == GLFW.KEY_LEFT || key == GLFW.KEY_A)
		right = normalized(cross(viewer.dir, viewer.vup))
		viewer.pos = viewer.pos - right * viewer.walk_speed
		redisplay(viewer)
		return
	end

	if (key == GLFW.KEY_RIGHT || key == GLFW.KEY_D)
		right = normalized(cross(viewer.dir, viewer.vup))
		viewer.pos = viewer.pos + right * viewer.walk_speed
		redisplay(viewer)
		return
	end



	if (key == GLFW.KEY_L)
		viewer.show_lines = !viewer.show_lines
		redisplay(viewer)
		return
	end


end

# ///////////////////////////////////////////////////////////////////////////////////////////////////////
# Reference to GKS ISO Graphics Standard (ISO/IEC 7942-4:1998)
export GLText

vectors = [
	([0.0, 0.0], [[1]]),
	([1.75 1.75 2.0 2.0 1.5 1.5; 1.75 5.5 0.5 1.0 1.0 0.5], [[1, 2], [3, 4], [4, 5], [5, 6], [6, 3]]),
	([1.0 2.0 2.0 2.0 1.5 1.5 2.0 3.0 3.0 3.0 2.5 2.5; 4.0 5.0 5.5 6.0 6.0 5.5 4.0 5.0 5.5 6.0 6.0 5.5], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 3], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12], [12, 9]]),
	([1.0 3.0 1.0 3.0 1.25 1.75 2.25 2.75; 2.5 2.5 3.5 3.5 1.75 4.0 1.75 4.0], [[1, 2], [3, 4], [5, 6], [7, 8]]),
	([0.0 1.0 3.0 4.0 4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0 2.0 2.0; 1.0 0.0 0.0 1.0 2.0 3.0 3.0 4.0 5.0 6.0 6.0 5.0 -0.5 6.5], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12], [13, 14]]),
	([2.5 2.5 2.0 2.0 2.5 2.5 2.0 2.0 1.5 3.5; 0.0 0.5 0.5 0.0 5.5 6.0 6.0 5.5 5.5 11.5], [[1, 2], [2, 3], [3, 4], [4, 1], [5, 6], [6, 7], [7, 8], [8, 5], [9, 10]]),
	([4.0 3.0 1.0 0.0 0.0 1.0 2.0 3.0 3.0 2.0 1.0 1.0 3.0 4.0; 1.0 0.0 0.0 1.0 2.0 3.0 3.0 4.0 5.0 6.0 5.0 4.0 2.0 2.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12], [12, 13], [13, 14]]),
	([1.0 2.0 2.0 2.0 1.5 1.5; 4.0 5.0 5.5 6.0 6.0 5.5], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 3]]),
	([2.0 1.0 0.5 1.0 2.0; 0.0 1.0 3.0 5.0 6.0], [[1, 2], [2, 3], [3, 4], [4, 5]]),
	([2.0 3.0 3.5 3.0 2.0; 0.0 1.0 3.0 5.0 6.0], [[1, 2], [2, 3], [3, 4], [4, 5]]),
	([1.0 3.0 2.0 2.0 1.0 3.0 1.0 3.0; 3.0 3.0 2.0 4.0 2.0 4.0 4.0 2.0], [[1, 2], [3, 4], [5, 6], [7, 8]]),
	([1.0 3.0 2.0 2.0; 3.0 3.0 2.0 4.0], [[1, 2], [3, 4]]),
	([1.0 2.0 2.0 2.0 1.5 1.5; -1.0 0.0 0.5 1.0 1.0 0.5], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 3]]),
	([1.0 3.0; 3.0 3.0], [[1, 2]]),
	([2.0 2.0 1.5 1.5 2.0; 0.0 0.5 0.5 0.0 0.0], [[1, 2], [2, 3], [3, 4], [4, 5]]),
	([1.0 3.0; 0.0 6.0], [[1, 2]]),
	([4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0; 1.0 0.0 0.0 1.0 5.0 6.0 6.0 5.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 1], [4, 8]]),
	([0.0 2.0 2.0 0.0 4.0; 4.0 6.0 0.0 0.0 0.0], [[1, 2], [2, 3], [4, 5]]),
	([0.0 0.0 1.0 3.0 4.0 4.0 0.0 4.0; 4.0 5.0 6.0 6.0 5.0 4.0 0.0 0.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8]]),
	([0.0 4.0 2.0 4.0 4.0 3.0 1.0 0.0 0.0; 6.0 6.0 4.0 2.0 1.0 0.0 0.0 1.0 2.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9]]),
	([4.0 0.0 4.0 4.0; 1.0 1.0 6.0 0.0], [[1, 2], [2, 3], [3, 4]]),
	([4.0 0.0 0.0 3.0 4.0 4.0 3.0 1.0 0.0 0.0; 6.0 6.0 4.0 4.0 3.0 1.0 0.0 0.0 1.0 2.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10]]),
	([4.0 1.0 0.0 0.0 1.0 3.0 4.0 4.0 3.0 1.0 0.0; 6.0 6.0 5.0 1.0 0.0 0.0 1.0 3.0 4.0 4.0 3.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11]]),
	([0.0 0.0 4.0 0.0; 5.0 6.0 6.0 0.0], [[1, 2], [2, 3], [3, 4]]),
	([1.0 3.0 4.0 4.0 3.0 1.0 0.0 0.0 4.0 4.0 3.0 1.0 0.0 0.0; 0.0 0.0 1.0 2.0 3.0 3.0 2.0 1.0 4.0 5.0 6.0 6.0 5.0 4.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 1], [6, 5], [5, 9], [9, 10], [10, 11], [11, 12], [12, 13], [13, 14], [14, 6]]),
	([0.0 3.0 4.0 4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0; 0.0 0.0 1.0 5.0 6.0 6.0 5.0 3.0 2.0 2.0 3.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11]]),
	([2.0 2.0 1.5 1.5 2.0 2.0 1.5 1.5; 1.0 1.5 1.5 1.0 3.0 3.5 3.5 3.0], [[1, 2], [2, 3], [3, 4], [4, 1], [5, 6], [6, 7], [7, 8], [8, 5]]),
	([2.0 2.0 1.5 1.5 1.0 2.0 2.0 2.0 1.5 1.5; 3.0 3.5 3.5 3.0 -0.5 0.5 1.0 1.5 1.5 1.0], [[1, 2], [2, 3], [3, 4], [4, 1], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 7]]),
	([3.0 0.0 3.0; 6.0 3.0 0.0], [[1, 2], [2, 3]]),
	([1.0 3.0 1.0 3.0; 2.5 2.5 3.5 3.5], [[1, 2], [3, 4]]),
	([1.0 4.0 1.0; 6.0 3.0 0.0], [[1, 2], [2, 3]]),
	([2.0 2.0 1.5 1.5 1.75 1.75 3.0 3.0 2.0 1.0 0.0 0.0; 0.0 0.5 0.5 0.0 1.0 2.75 4.0 5.0 6.0 6.0 5.0 4.0], [[1, 2], [2, 3], [3, 4], [4, 1], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12]]),
	([4.0 1.0 0.0 0.0 1.0 3.0 4.0 4.0 2.0 1.0 2.0 3.0 2.0; 0.0 0.0 1.0 3.0 4.0 4.0 3.0 1.0 1.0 2.0 3.0 2.0 1.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12], [12, 13]]),
	([0.0 0.0 1.0 3.0 4.0 4.0 0.0 4.0; 0.0 5.0 6.0 6.0 5.0 0.0 2.0 2.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [7, 8]]),
	([0.0 0.0 3.0 4.0 4.0 3.0 4.0 4.0 3.0 0.0 0.0 3.0; 0.0 6.0 6.0 5.0 4.0 3.0 2.0 1.0 0.0 0.0 3.0 3.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [11, 12]]),
	([4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0; 1.0 0.0 0.0 1.0 5.0 6.0 6.0 5.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8]]),
	([0.0 0.0 3.0 4.0 4.0 3.0 0.0; 0.0 6.0 6.0 5.0 1.0 0.0 0.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]]),
	([4.0 0.0 0.0 4.0 0.0 3.0; 0.0 0.0 6.0 6.0 3.0 3.0], [[1, 2], [2, 3], [3, 4], [5, 6]]),
	([0.0 0.0 4.0 0.0 3.0; 0.0 6.0 6.0 3.0 3.0], [[1, 2], [2, 3], [4, 5]]),
	([2.0 4.0 4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0; 3.0 3.0 1.0 0.0 0.0 1.0 5.0 6.0 6.0 5.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10]]),
	([0.0 0.0 4.0 4.0 0.0 4.0; 0.0 6.0 6.0 0.0 3.0 3.0], [[1, 2], [3, 4], [5, 6]]),
	([2.0 2.0 1.0 3.0 1.0 3.0; 0.0 6.0 0.0 0.0 6.0 6.0], [[1, 2], [3, 4], [5, 6]]),
	([0.0 1.0 2.0 3.0 3.0 2.0 4.0; 1.0 0.0 0.0 1.0 6.0 6.0 6.0], [[1, 2], [2, 3], [3, 4], [4, 5], [6, 7]]),
	([4.0 0.0 4.0 0.0 0.0; 6.0 3.0 0.0 0.0 6.0], [[1, 2], [2, 3], [4, 5]]),
	([4.0 0.0 0.0; 0.0 0.0 6.0], [[1, 2], [2, 3]]),
	([0.0 0.0 2.0 4.0 4.0; 0.0 6.0 4.0 6.0 0.0], [[1, 2], [2, 3], [3, 4], [4, 5]]),
	([0.0 0.0 4.0 4.0 4.0; 0.0 6.0 2.0 0.0 6.0], [[1, 2], [2, 3], [4, 5]]),
	([4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0 4.0; 1.0 0.0 0.0 1.0 5.0 6.0 6.0 5.0 1.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9]]),
	([0.0 0.0 3.0 4.0 4.0 3.0 0.0; 0.0 6.0 6.0 5.0 3.0 2.0 2.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]]),
	([4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0 4.0 3.0 4.0; 1.0 0.0 0.0 1.0 5.0 6.0 6.0 5.0 1.0 1.0 0.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [10, 11]]),
	([0.0 0.0 3.0 4.0 4.0 3.0 0.0 3.0 4.0; 0.0 6.0 6.0 5.0 3.0 2.0 2.0 2.0 0.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [8, 9]]),
	([0.0 1.0 3.0 4.0 4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0; 1.0 0.0 0.0 1.0 2.0 3.0 3.0 4.0 5.0 6.0 6.0 5.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12]]),
	([2.0 2.0 0.0 4.0; 0.0 6.0 6.0 6.0], [[1, 2], [3, 4]]),
	([0.0 0.0 1.0 3.0 4.0 4.0; 6.0 1.0 0.0 0.0 1.0 6.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6]]),
	([0.0 2.0 4.0; 6.0 0.0 6.0], [[1, 2], [2, 3]]),
	([0.0 0.0 1.0 2.0 3.0 4.0 4.0; 6.0 3.0 0.0 3.0 0.0 3.0 6.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]]),
	([0.0 4.0 0.0 4.0; 0.0 6.0 6.0 0.0], [[1, 2], [3, 4]]),
	([0.0 2.0 4.0 2.0 2.0; 6.0 2.0 6.0 2.0 0.0], [[1, 2], [2, 3], [4, 5]]),
	([0.0 4.0 0.0 4.0; 6.0 6.0 0.0 0.0], [[1, 2], [2, 3], [3, 4]]),
	([2.0 1.0 1.0 2.0; 0.0 0.0 6.0 6.0], [[1, 2], [2, 3], [3, 4]]),
	([1.0 3.0; 6.0 0.0], [[1, 2]]),
	([2.0 3.0 3.0 2.0; 0.0 0.0 6.0 6.0], [[1, 2], [2, 3], [3, 4]]),
	([1.0 2.0 3.0; 5.0 6.0 5.0], [[1, 2], [2, 3]]),
	([1.0 4.0; 0.0 0.0], [[1, 2]]),
	([2.0 2.0 3.0 2.5 2.5 2.0; 4.5 5.0 6.0 4.0 4.5 4.0], [[1, 2], [2, 3], [4, 5], [5, 1], [1, 6], [6, 4]]),
	([4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0 4.0 4.0; 1.0 0.0 0.0 1.0 2.0 3.0 3.0 2.0 0.0 3.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 1], [9, 10]]),
	([4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0 0.0 0.0 1.0; 1.0 0.0 0.0 1.0 2.0 3.0 3.0 2.0 0.0 5.0 5.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 1], [9, 10], [10, 11]]),
	([4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0; 1.0 0.0 0.0 1.0 2.0 3.0 3.0 2.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8]]),
	([4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0 4.0 4.0 3.0; 1.0 0.0 0.0 1.0 2.0 3.0 3.0 2.0 0.0 5.0 5.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 1], [9, 10], [10, 11]]),
	([4.0 1.0 0.0 0.0 1.0 3.0 4.0 4.0 0.0; 0.0 0.0 1.0 2.0 3.0 3.0 2.0 1.0 1.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9]]),
	([4.0 4.0 3.0 2.0 1.0 1.0 0.0 2.0; 3.0 4.0 5.0 5.0 4.0 0.0 1.0 1.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [7, 8]]),
	([4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0 4.0 3.0 1.0 0.0; 1.0 0.0 0.0 1.0 2.0 3.0 3.0 2.0 0.0 -1.0 -1.0 0.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 1], [1, 9], [9, 10], [10, 11], [11, 12]]),
	([4.0 4.0 3.0 1.0 0.0 0.0 0.0 1.0; 0.0 2.0 3.0 3.0 2.0 0.0 5.0 5.0], [[1, 2], [2, 3], [3, 4], [4, 5], [6, 7], [7, 8]]),
	([1.0 3.0 1.0 3.0 2.0 2.0 2.25 2.25 1.75 1.75; 0.0 0.0 3.0 3.0 0.0 3.0 3.75 4.25 4.25 3.75], [[1, 2], [3, 4], [5, 6], [7, 8], [8, 9], [9, 10], [10, 7]]),
	([1.0 3.0 2.0 2.0 1.0 0.0 2.25 2.25 1.75 1.75; 3.0 3.0 3.0 0.0 -1.0 0.0 3.75 4.25 4.25 3.75], [[1, 2], [3, 4], [4, 5], [5, 6], [7, 8], [8, 9], [9, 10], [10, 7]]),
	([0.0 1.0 1.0 0.0 4.0 2.0 1.0 3.0 4.0; 0.0 0.0 3.0 3.0 0.0 0.0 1.0 3.0 3.0], [[1, 2], [2, 3], [3, 4], [5, 6], [6, 7], [7, 8], [8, 9]]),
	([2.0 2.0 1.0 1.0 3.0; 0.0 5.0 5.0 0.0 0.0], [[1, 2], [2, 3], [4, 5]]),
	([4.0 4.0 2.0 2.0 2.0 0.0 0.0 0.0; 0.0 3.0 2.0 0.0 3.0 2.0 0.0 3.0], [[1, 2], [2, 3], [4, 5], [5, 6], [7, 8]]),
	([3.0 3.0 1.0 1.0 1.0; 0.0 3.0 2.0 0.0 3.0], [[1, 2], [2, 3], [4, 5]]),
	([4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0 4.0; 1.0 0.0 0.0 1.0 2.0 3.0 3.0 2.0 1.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9]]),
	([4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0 0.0 0.0; 1.0 0.0 0.0 1.0 2.0 3.0 3.0 2.0 3.0 -1.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 1], [9, 10]]),
	([4.0 3.0 1.0 0.0 0.0 1.0 3.0 4.0 4.0 4.0; 1.0 0.0 0.0 1.0 2.0 3.0 3.0 2.0 3.0 -1.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 1], [9, 10]]),
	([0.0 2.0 1.0 1.0 1.0 2.0 3.0 4.0; 0.0 0.0 0.0 3.0 2.0 3.0 3.0 2.0], [[1, 2], [3, 4], [5, 6], [6, 7], [7, 8]]),
	([0.0 4.0 3.0 1.0 0.0 1.0 3.0 4.0; 0.0 0.0 1.0 1.0 2.0 3.0 3.0 2.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8]]),
	([1.0 3.0 2.0 2.0 2.0 3.0; 0.0 0.0 0.0 5.0 4.0 4.0], [[1, 2], [3, 4], [5, 6]]),
	([0.0 1.0 1.0 2.0 3.0 4.0 4.0; 3.0 3.0 1.0 0.0 0.0 1.0 3.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]]),
	([0.0 1.0 3.0 4.0; 3.0 0.0 3.0 3.0], [[1, 2], [2, 3], [3, 4]]),
	([0.0 0.0 1.0 2.0 3.0 4.0 4.0; 3.0 2.0 0.0 2.0 0.0 2.0 3.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]]),
	([0.0 1.0 4.0 1.0 4.0; 3.0 3.0 0.0 0.0 3.0], [[1, 2], [2, 3], [4, 5]]),
	([0.0 1.0 2.5 0.0 1.0 4.0; 3.0 3.0 1.5 0.0 0.0 3.0], [[1, 2], [2, 3], [4, 5], [5, 6]]),
	([0.0 0.0 3.0 0.0 3.0 4.0; 2.0 3.0 3.0 0.0 0.0 1.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6]]),
	([2.5 2.0 2.0 1.5 2.0 2.0 2.5; 6.5 6.0 3.5 3.0 2.5 0.0 -0.5], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]]),
	([2.0 2.0; 0.0 5.0], [[1, 2]]),
	([1.5 2.0 2.0 2.5 2.0 2.0 1.5; 6.5 6.0 3.5 3.0 2.5 0.0 -0.5], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]]),
	([1.0 1.75 2.75 3.5; 5.0 5.5 5.0 5.5], [[1, 2], [2, 3], [3, 4]])]

# seems coordinates are in the range 4x6, is it right?


# ////////////////////////////////////////////////////////////////////////////
function GLText(mystring; center=Point3d(0.0, 0.0, 0.0), align="center", fontsize=0.05, color=Point4d(0, 0, 0, 1))

	ret = Vector{GLBatch}()

	if mystring == ""
		return ret
	end

	sx = DEFAULT_TEXT_SCALING[1] * fontsize
	sy = DEFAULT_TEXT_SCALING[2] * fontsize

	vertices = Vector{Float32}()
	colors = Vector{Float32}()
	tx, ty = 0.0, 0.0
	x1, y1, z1 = +Inf, +Inf, +Inf
	x2, y2, z2 = -Inf, -Inf, -Inf
	for ch in mystring
		coords, indices = vectors[1+Int(ch)-32]
		for (A, B) in indices

			p0 = Point3d(coords[1, A] * sx + tx, coords[2, A] * sy + ty, 0.0)
			p1 = Point3d(coords[1, B] * sx + tx, coords[2, B] * sy + ty, 0.0)
			append!(vertices, p0)
			append!(colors, color)
			append!(vertices, p1)
			append!(colors, color)
			x1 = min(x1, p0[1], p1[1])
			y1 = min(y1, p0[2], p1[2])
			z1 = min(z1, p0[3], p1[3])
			x2 = max(x2, p0[1], p1[1])
			y2 = max(y2, p0[2], p1[2])
			z2 = max(z2, p0[3], p1[3])
		end
		tx += fontsize
		ty += 0
	end

	# println(x1," ",x2," ",y1," ", y2," ",z1,z2)
	vt = center

	# assuming center alignmet
	if true
		for I in 1:3:length(vertices)
			vertices[I+0] = vertices[I+0] + center[1] - (x2 - x1) / 2.0
			vertices[I+1] = vertices[I+1] + center[2] - (y2 - y1) / 2.0
			vertices[I+2] = vertices[I+2] + center[3] - (z2 - z1) / 2.0
		end
	end

	batch = GLBatch(LINES)
	batch.vertices = GLVertexBuffer(vertices)
	batch.colors = GLVertexBuffer(colors)

	push!(ret, batch)
	return ret
end

include("./defaults.jl")

# ///////////////////////////////////////////////////////////////////
function GLText(
	V::Vector{Vector{Float64}};
	EV=nothing,
	FV=nothing,
	V_color=WHITE,
	EV_color=BLUE,
	FV_color=GREEN,
	fontsize=0.05)

	batches = Vector{GLBatch}()

	if V_color[4] > 0
		for I in 1:length(V)
			append!(batches, GLText(string(I), center=ComputeCentroid([V[it] for it in [I]]), fontsize=fontsize, color=V_color))
		end
	end

	if EV != nothing && EV_color[4] > 0
		for I in 1:length(EV)
			append!(batches, GLText(string(I), center=ComputeCentroid([V[it] for it in EV[I]]), fontsize=fontsize, color=EV_color))
		end
	end

	if FV != nothing && FV_color[4] > 0
		for I in 1:length(FV)
			append!(batches, GLText(string(I), center=ComputeCentroid([V[it] for it in FV[I]]), fontsize=fontsize, color=FV_color))
		end
	end

	return batches

end

# ///////////////////////////////////////////////////////////////

# ///////////////////////////////////////////////////////////////////////////
# scrgiorgio: not sure this class is needed
mutable struct GLColor

	choice::Int   # to select between per_body, per_cell, and per_vertex

	color_per_body::Vector{Point4d}
	color_per_cell::Vector{Point4d}
	color_per_vertex::Vector{Point4d}

	colors::GLVertexBuffer

	# constructor
	GLColor() = GLColor(1, COLORS[1], Point4d[], Point4d[], Point4d[], Point4d[])

	# constructor
	GLColor(choice) = GLColor(choice, COLORS[1], Point4d[], Point4d[], Point4d[], Point4d[])

	# constructor
	function GLColor(choice, color_per_body, color_per_cell, color_per_vertex, colors)
		ret = new(choice, color_per_body, color_per_cell, color_per_vertex, colors)
		return ret
	end

end




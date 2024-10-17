using LinearAlgebra
using StaticArrays

# /////////////////////////////////////////////////////////////////////
import Base: *
import Base.:-
import Base.:+

const ConcretePointNd  = Vector{Float64}
const ConcretePointsNd = Vector{ConcretePointNd}

const AbstractPointNd  = Vector{<:Number}
const AbstractPointsNd = Vector{<:AbstractPointNd}

const PointNd  = ConcretePointNd
const PointsNd = ConcretePointsNd

function to_concrete(value::ConcretePointNd)::ConcretePointNd
	return value
end

function to_concrete(value::ConcretePointsNd)::ConcretePointsNd
	return value
end

function to_concrete(value::AbstractPointNd)::ConcretePointNd
	return PointNd([convert(Float64,it) for it in value])
end

function to_concrete(value::AbstractPointsNd)::ConcretePointsNd
	return PointsNd([to_concrete(it) for it in value])
end

export AbstractPointNd
export AbstractPointsNd
export PointNd
export PointsNd

# //////////////////////////////////////////////////////////////////////////////
const Point3d = MVector{3,Float64}
export Point3d

Point3d() = Point3d(0.0, 0.0, 0.0)


function norm(p::Point3d)
	return sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3])
end
export norm

function normalized(p::Point3d)
	len = norm(p)
	return Point3d(p[1] / len, p[2] / len, p[3] / len)
end
export normalized

function +(a::Point3d, b::Point3d)
	return Point3d(a[1] + b[1], a[2] + b[2], a[3] + b[3])
end

function -(a::Point3d, b::Point3d)
	return Point3d(a[1] - b[1], a[2] - b[2], a[3] - b[3])
end

function *(a::Point3d, vs::Float64)
	return Point3d(a[1] * vs, a[2] * vs, a[3] * vs)
end

# ///////////////////////////////////////////////////////////////////////
function computeNormal(p0::Point3d, p1::Point3d, p2::Point3d)
	return normalized(cross(p1 - p0, p2 - p0))
end
export computeNormal


# //////////////////////////////////////////////////////////////////////////////
const Point4d = MVector{4,Float64}
export Point4d

Point4d() = Point4d(0.0, 0.0, 0.0, 0.0)

function norm(p::Point4d)
	return sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + p[4] * p[4])
end

function normalized(p::Point4d)
	len = norm(p)
	return Point4d(p[1] / len, p[2] / len, p[3] / len, p[4] / len)
end

function dropW(p::Point4d)
	return Point4d(p[1] / p[4], p[2] / p[4], p[3] / p[4])
end

function +(a::Point4d, b::Point4d)
	return Point4d(a[1] + b[1], a[2] + b[2], a[3] + b[3], a[4] + b[4])
end

function -(a::Point4d, b::Point4d)
	return Point4d(a[1] - b[1], a[2] - b[2], a[3] - b[3], a[4] - b[4])
end

function *(a::Point4d, vs::Float64)
	return Point4d(a[1] * vs, a[2] * vs, a[3] * vs, a[4] * vs)
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
export Box3d

Base.:(==)(a::Box3d, b::Box3d) = a.p1 == b.p1 && a.p2 == b.p2

function invalidBox()::Box3d
	m, M = typemin(Float64), typemax(Float64)
	return Box3d(Point3d(M, M, M), Point3d(m, m, m))
end
export invalidBox

function addPoint(box::Box3d, point::Point3d)
	for i in 1:3
		box.p1[i] = min(box.p1[i], point[i])
		box.p2[i] = max(box.p2[i], point[i])
	end
	return box
end

function getPoints(box::Box3d)::Vector{Point3d}
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
export getPoints

function center(box::Box3d)::Point3d
	return Point3d(
		0.5 * box.p1[1] + 0.5 * box.p2[1],
		0.5 * box.p1[2] + 0.5 * box.p2[2],
		0.5 * box.p1[3] + 0.5 * box.p2[3])
end
export center

const Matrix3d = MMatrix{3,3,Float64}
export Matrix3d

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
export flatten

# ////////////////////////////////////////////////////////////////////
const Matrix4d = MMatrix{4,4,Float64}
export Matrix3d

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


function Matrix4d(v::PointNd)
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
export translateMatrix

function scaleMatrix(vs::Point3d)
	return Matrix4d(
		vs[1], 0.0, 0.0, 0.0,
		0.0, vs[2], 0.0, 0.0,
		0.0, 0.0, vs[3], 0.0,
		0.0, 0.0, 0.0, 1.0)
end
export scaleMatrix

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
export lookAtMatrix


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
export perspectiveMatrix


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
export orthoMatrix

function getLookAt(T::Matrix4d)
	T = inv(T)
	pos = (Point3d(T[1, 4], T[2, 4], T[3, 4]))
	dir = normalized(Point3d(-T[1, 3], -T[2, 3], -T[3, 3]))
	vup = normalized(Point3d(T[1, 2], T[2, 2], T[3, 2]))
	pos, dir, vup
end
export getLookAt

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
export Quaternion

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
export convertToQuaternion

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

export convertToMatrix


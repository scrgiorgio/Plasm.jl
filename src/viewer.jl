Properties = Dict{String,Any}
export Properties

# ///////////////////////////////////////////////////////////////////////
function split(v::AbstractVector, size::Int)::AbstractVector
	N=length(v)
	@assert(Int(N/size)*size==N)
	ret=[]
	for C in 1:size:N
		push!(ret, copy(collect(v[C:C+size-1])))
	end
	return ret
end
export split

# ///////////////////////////////////////////////////////////////////////
function transform_points(T::Matrix4d, points::Vector{Float32})::Vector{Float32}
	ret=Vector{Float32}()
	for (x,y,z) in split(points,3)
		x,y,z,w=[
			T[1,1]*x + T[1,2]*y + T[1,3]*z +  T[1,4]*1.0,
			T[2,1]*x + T[2,2]*y + T[2,3]*z +  T[2,4]*1.0,
			T[3,1]*x + T[3,2]*y + T[3,3]*z +  T[3,4]*1.0,
			T[4,1]*x + T[4,2]*y + T[4,3]*z +  T[4,4]*1.0
		]
		append!(ret,[Float32(x/w),Float32(y/w),Float32(z/w)]...)
	end
	return ret
end
export transform_points

# ///////////////////////////////////////////////////////////////////////
function triangles_to_lines(points::Vector{Float32})::Vector{Float32}
	ret=Vector{Float32}()
	for (x0,y0,z0,x1,y1,z1,x2,y2,z2) in split(points,9)
		p0=Point3d(x0,y0,z0)
		p1=Point3d(x1,y1,z1)
		p2=Point3d(x2,y2,z2)
		append!(ret,p0...);append!(ret,p1...)
		append!(ret,p1...);append!(ret,p2...)
		append!(ret,p2...);append!(ret,p0...)
	end
	return ret
end
export triangles_to_lines

# ///////////////////////////////////////////////////////////////////////
function compute_triangles_normals(points::Vector{Float32})::Vector{Float32}
	ret=Vector{Float32}()
	for (x0,y0,z0, x1,y1,z1, x2,y2,z2) in split(points,9)
		p0=Point3d(x0,y0,z0)
		p1=Point3d(x1,y1,z1)
		p2=Point3d(x2,y2,z2)
		normal = computeNormal(p0, p1, p2)
		append!(ret,normal...)
		append!(ret,normal...)
		append!(ret,normal...)
	end
	return ret
end
export compute_triangles_normals

include("./viewer.glfw.jl")
# include("./viewer.meshcat.jl")


# ////////////////////////////////////////////////////////////////////////
function render_cuboid(viewer::Viewer, box::Box3d)
	box_points::Vector{Point3d} = getPoints(box)
	faces = [[1, 2, 3, 4], [4, 3, 7, 8], [8, 7, 6, 5], [5, 6, 2, 1], [6, 7, 3, 2], [8, 5, 1, 4]]
	points  = Vector{Float32}()
	for face in faces
		p3, p2, p1, p0 = box_points[face[1]], box_points[face[2]], box_points[face[3]], box_points[face[4]] # reverse order
		append!(points, p0);append!(points, p1);append!(points, p2)
		append!(points, p0);append!(points, p2);append!(points, p3)
	end
	render_triangles(viewer, points)
	render_lines(viewer, triangles_to_lines(points), line_width=1)
end
export render_cuboid 


# ////////////////////////////////////////////////////////////////////////
function render_axis(viewer::Viewer, p0::Point3d, p1::Point3d)
	vertices = Vector{Float32}()
	colors   = Vector{Float32}()
	R = Point4d(1, 0, 0, 1)
	G = Point4d(0, 1, 0, 1)
	B = Point4d(0, 0, 1, 1)
	append!(vertices, p0);append!(vertices, Point3d(p1[1], p0[2], p0[3]));append!(colors, R);append!(colors, R)
	append!(vertices, p0);append!(vertices, Point3d(p0[1], p1[2], p0[3]));append!(colors, G);append!(colors, G)
	append!(vertices, p0);append!(vertices, Point3d(p0[1], p0[2], p1[3]));append!(colors, B);append!(colors, B)
	render_lines(viewer, vertices, colors=colors, line_width=2)
end

export render_axis


# ///////////////////////////////////////////////////////////////////////////////////////////////////////
# Reference to GKS ISO Graphics Standard (ISO/IEC 7942-4:1998)
gltext_char_vectors = [
	
	([0.0, 0.0], [[1,1]]), # space
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
function render_text(viewer::Viewer, s::String; center=Point3d(0.0, 0.0, 0.0), align="center", fontsize=0.05, color=Point4d(0, 0, 0, 1))

	if s == ""
		return 
	end

	sx = DEFAULT_TEXT_SCALING[1] * fontsize
	sy = DEFAULT_TEXT_SCALING[2] * fontsize

	vertices = Vector{Float32}()
	colors = Vector{Float32}()
	tx, ty = 0.0, 0.0
	x1, y1, z1 = +Inf, +Inf, +Inf
	x2, y2, z2 = -Inf, -Inf, -Inf
	for ch in s
		coords, indices = gltext_char_vectors[1+Int(ch)-32]
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
	begin
		for I in 1:3:length(vertices)
			vertices[I+0] = vertices[I+0] + center[1] - (x2 - x1) / 2.0
			vertices[I+1] = vertices[I+1] + center[2] - (y2 - y1) / 2.0
			vertices[I+2] = vertices[I+2] + center[3] - (z2 - z1) / 2.0
		end
	end

	render_lines(viewer, vertices, colors=colors, line_width=1)
end

export render_text
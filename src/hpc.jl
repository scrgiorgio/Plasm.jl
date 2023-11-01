using Test,LinearAlgebra,PyCall

export ComputeTriangleNormal,GoodTetOrientation,
	BoxNd,MatrixNd,Hpc,BuildMkPol,
	toList,valid,fuzzyEqual,dim,size,center,addPoint,addPoints,addBox,isIdentity,transpose,invert,dim,embed,adjoin,transformPoint,translate,scale,rotate,box,
	MkPol,Struct,Cube,Simplex,Join,Quote,Transform,Translate,Scale,Rotate,Power,UkPol,GetBatches,MapFn,
	ToSimplicialForm,ToBoundaryForm,ToLARForm,
	View

import Base.:(==)
import Base.:*
import Base.size
import Base.transpose

# /////////////////////////////////////////////////////////////
function ComputeTriangleNormal(p0::Vector{Float64}, p1::Vector{Float64}, p2::Vector{Float64})
	 p0 = vcat(p0, zeros(3 - length(p0)))
	 p1 = vcat(p1, zeros(3 - length(p1)))
	 p2 = vcat(p2, zeros(3 - length(p2)))
	 a = [p1[i] - p0[i] for i in 1:3]
	 b = [p2[i] - p0[i] for i in 1:3]
	 ret = LinearAlgebra.cross(a, b)
	 N = LinearAlgebra.norm(ret)
	 if N == 0.0
		  N = 1.0
	 end
	 return [ret[i] / N for i in 1:3]
end

# /////////////////////////////////////////////////////////////
function GoodTetOrientation(v0::Vector{Float64}, v1::Vector{Float64}, v2::Vector{Float64}, v3::Vector{Float64})
	 v0 = vcat(v0, zeros(3 - length(v0)))
	 v1 = vcat(v1, zeros(3 - length(v1)))
	 v2 = vcat(v2, zeros(3 - length(v2)))
	 v3 = vcat(v3, zeros(3 - length(v3)))
	 a = [v3[i] - v1[i] for i in 1:3]
	 b = [v2[i] - v1[i] for i in 1:3]
	 c = [v0[i] - v1[i] for i in 1:3]
	 n = LinearAlgebra.cross(a, b)
	 return dot(n, c) > 0
end

# /////////////////////////////////////////////////////////////
mutable struct BoxNd
	 p1::Vector{Float64}
	 p2::Vector{Float64}

	function BoxNd(dim::Int)
		p1 = [+floatmax(Float64) for I in 1:dim]
		p2 = [-floatmax(Float64) for I in 1:dim]	
		new(p1, p2)
	end
	
	function BoxNd(p1::Vector{Float64},p2::Vector{Float64})
		@assert length(p1) == length(p2)
		p1 = copy(p1)
		p2 = copy(p2)
		new(p1, p2)
	end

end


toList(self::BoxNd) = [copy(self.p1), copy(self.p2)]

function valid(self::BoxNd)
	 for i in 1:dim(self)
		  if self.p1[i] > self.p2[i]
				return false
		  end
	 end
	 return true
end

==(box1::BoxNd, box2::BoxNd) = isa(box1, typeof(box2)) && box1.p1 == box2.p1 && box1.p2 == box2.p2

function fuzzyEqual(box1::BoxNd, box2::BoxNd, Epsilon=1e-4)
	 if !(isa(box2, typeof(box1)) && dim(box1) == dim(box2))
		  return false
	 end
	 p1 = [abs(a - b) <= Epsilon for (a, b) in zip(box1.p1, box2.p1)]
	 p2 = [abs(a - b) <= Epsilon for (a, b) in zip(box1.p2, box2.p2)]
	 return !(false in p1) && !(false in p2)
end

Base.show(io::IO, self::BoxNd) = print(io, "BoxNd(", repr(self.p1), ", ", repr(self.p2), ")")

dim(self::BoxNd) = length(self.p1)

function size(self::BoxNd)
	 return [To - From for (From, To) in zip(self.p1, self.p2)]
end

function center(self::BoxNd)
	 return [0.5 * (From + To) for (From, To) in zip(self.p1, self.p2)]
end

function addPoint(self::BoxNd, point::Vector{Float64})
	 for i in 1:dim(self)
		self.p1[i] = min(self.p1[i], point[i])
		self.p2[i] = max(self.p2[i], point[i])
	 end
	 return self
end

function addPoints(self::BoxNd, points)
	 for point in points
		  addPoint(self, point)
	 end
	 return self
end

function addBox(box1::BoxNd, box2::BoxNd)
	 return addPoint(box1, box2.p1).addPoint(box2.p2)
end

# /////////////////////////////////////////////////////////////
mutable struct MatrixNd
	T::Matrix{Float64}

	function MatrixNd(dim::Int=0)
		T = Matrix{Float64}(I, dim, dim)
		new(T)
	end

	function MatrixNd(other::MatrixNd)
		new(copy(other.T))
	end	
	
	function MatrixNd(T::Matrix{Float64})
		new(copy(T))
	end

	function MatrixNd(arg::Vector{Vector{Float64}})
		T = reduce(vcat, arg')
		new(T)
	end

	function MatrixNd(arg::Vector{Vector{Int64}})
		T = reduce(vcat, arg')
		new(T)
	end	

end


Base.getindex(self::MatrixNd, args...) = getindex(self.T, args...)
Base.setindex!(self::MatrixNd, args...) = setindex!(self.T, args...)

==(matrix1::MatrixNd, matrix2::MatrixNd) = isa(matrix1, typeof(matrix2)) && matrix1.T == matrix2.T

Base.show(io::IO, self::MatrixNd) = print(io, "MatrixNd(", repr(toList(self)), ")")

function isIdentity(self::MatrixNd)
	 return self.T == I
end

toList(self::MatrixNd) = [self.T[R,:] for R in 1:size(self.T,1)]

function transpose(self::MatrixNd)
	 return MatrixNd(Matrix{Float64}(transpose(self.T)))
end

function invert(self::MatrixNd)
	 return MatrixNd(inv(self.T))
end

dim(self::MatrixNd) = size(self.T, 1)

function embed(self::MatrixNd, target_dim)
	current_dim=dim(self)
	 if target_dim <= current_dim
		  return self
	 end
	 ret = MatrixNd(target_dim)
	 ret.T[1:current_dim, 1:current_dim] = self.T
	 return ret
end

function adjoin(matrix1::MatrixNd, matrix2::MatrixNd)
	M, N = dim(matrix1), dim(matrix2)
	T=M + N - 1
	ret = MatrixNd(T)
	ret[2:M, 2:M] = matrix1[2:M, 2:M]
	for I in 2:M
		ret[I, 1] = matrix1[I, 1]
		ret[1, I] = matrix1[1, I]
	end
	ret[M+1:T, M+1:T] = matrix2[2:N, 2:N]
	for I in 2:N
		ret[I+M-1, 1] = matrix2[I, 1]
		ret[1, I+M-1] = matrix2[1, I]
	end
	return ret
end

function *(matrix1::MatrixNd, matrix2::MatrixNd)
	 return MatrixNd(matrix1.T * matrix2.T)
end

function transformPoint(self::MatrixNd, point::Vector{Float64})
	 point = self.T * [1.0; point; zeros(dim(self) - length(point) - 1)]
	 return [point[i, 1] / point[1, 1] for i in 2:dim(self)]
end

function translate(vt)
	 T = MatrixNd(length(vt) + 1)
	 for I in 2:dim(T)
		  T[I, 1] = vt[I-1]
	 end
	 return T
end

function scale(vs)
	 T = MatrixNd(length(vs) + 1)
	 for I in 2:dim(T)
		  T[I, I] = vs[I-1]
	 end
	 return T
end

function rotate(i, j, angle)
	i+=1
	j+=1
	T = MatrixNd(max(i, j))
	T[i, i] = +cos(angle)
	T[i, j] = -sin(angle)
	T[j, i] = +sin(angle)
	T[j, j] = +cos(angle)
	return T
end


# /////////////////////////////////////////////////////////////
mutable struct BuildMkPol

	 db::Dict{Vector{Float64}, Int}
	 points::Vector{Vector{Float64}}
	 hulls::Vector{Vector{Int}}
	 facets::Vector{Vector{Int}}

	# constructor
	function BuildMkPol()
		self=new(
			# db
			Dict{Vector{Float64}, Int}(),

			# points
			Vector{Vector{Float64}}(),

			# hulls
			Vector{Vector{Int}}(),

			# facets (for LAR)
			Vector{Vector{Int}}(),
		)
		return self
	end

	# constructor
	function BuildMkPol(points::Vector{Vector{Float64}}, hulls::Vector{Vector{Int}}=Vector{Vector{Int}}())

		self=BuildMkPol()

		if isempty(points)
			 return self
		end
		
		if isempty(hulls)
			 hulls = [collect(1:length(points))]
		end

		pdim = length(points[1])
		# special case in zero-dimension
		if pdim==0
			self.points=Vector{Vector{Float64}}()
			push!(self.points,[]) # an 'empty' point
			self.hulls=[[1]]
			return self
		end

		for hull in hulls

			 if isempty(hull)
				  continue
			 end

			 if pdim == 1
				  @assert length(hull) >= 2
				  box = BoxNd(1)
				  for idx in hull
						box = addPoint(box, points[idx])
				  end
				  addHull(self, Vector{Vector{Float64}}([box.p1, box.p2]))
			 else
				  hull_points = [points[idx] for idx in hull]
				  __spatial = pyimport_conda("scipy.spatial", "scipy") # the second argument is the conda package name
				  ConvexHull = __spatial.ConvexHull
				  try
						h = ConvexHull([points[idx] for idx in hull])
						hull_points = [Vector{Float64}(h.points[idx+1,:]) for idx in h.vertices]
				  catch
				  end
				  addHull(self, hull_points)
			 end
		end

		return self

  end	

end


function addPoint(self::BuildMkPol, p::Vector{Float64})::Int
	 idx = get(self.db, p, 0)
	 if idx>=1
		  return idx
	 else
		  idx = length(self.db) + 1
		  self.db[p] = idx
		  push!(self.points, p)
		  return idx
	 end
end


function addPoints(self::BuildMkPol, points::Vector{Vector{Float64}})::Dict{Int,Int}
	ret = Dict{Int,Int}()
	for P in 1:length(points)
		ret[P]=addPoint(self,points[P])
	end
	return ret
end

function addHull(self::BuildMkPol, points::Vector{Vector{Float64}})
	push!(self.hulls, [addPoint(self, p) for p in points])
end

dim(self::BuildMkPol) = isempty(self.points) ? 0 : length(self.points[1])

function box(self::BuildMkPol)
	 ret = BoxNd(dim(self))
	 if !isempty(self.points)
		  addPoints(ret, self.points)
	 end
	 return ret
end

Base.show(io::IO, self::BuildMkPol) = print(io, "BuildMkPol(points=", repr(self.points), ", hulls=", repr(self.hulls), ")")

# /////////////////////////////////////////////////////////////////////////////////
function ToSimplicialForm(self::BuildMkPol)

	 if isempty(self.points) || isempty(self.hulls)
		  return self
	 end
	 pdim = dim(self)
	 
	 if pdim <= 1
		  return self
	 end
	 
	 ret = BuildMkPol()
	 for hull in self.hulls
		  if length(hull) <= pdim + 1
				addHull(ret, [self.points[idx] for idx in hull])
		  else
				
				__spatial = pyimport_conda("scipy.spatial", "scipy") # the second argument is the conda package name
				ConvexHull = __spatial.ConvexHull
				Delaunay   = __spatial.Delaunay			
				try
					 d = Delaunay([self.points[idx] for idx in hull])
					 for simplex in [ d.simplices[R,:] for R in 1:size(d.simplices,1)]
							simplex_points = [Vector{Float64}(d.points[idx+1,:]) for idx in simplex]
						  addHull(ret, simplex_points)
					 end
				catch 
				end
		  end
	 end
	 FixOrientation!(ret)
	 return ret
end




# ////////////////////////////////////////////////////////////////////////////////
function FixOrientation!(self::BuildMkPol)
	 pdim = dim(self)
	 if pdim != 2 && pdim != 3
		  return
	 end
	 
	 if pdim == 2
		  fixed = Vector{Vector{Int}}()
		  for simplex in self.hulls
				if length(simplex) == 3
					 p0 = self.points[simplex[1]]
					 p1 = self.points[simplex[2]]
					 p2 = self.points[simplex[3]]
					 n = ComputeTriangleNormal(p0, p1, p2)
					 if n[3] < 0
						  simplex = [simplex[3], simplex[2], simplex[1]]
					 end
				end
				push!(fixed, simplex)
		  end
		  self.hulls = fixed
	 end
	 
	 if pdim == 3
		  fixed = Vector{Vector{Int}}()
		  for simplex in self.hulls
				if length(simplex) == 4
					 p0 = self.points[simplex[1]]
					 p1 = self.points[simplex[2]]
					 p2 = self.points[simplex[3]]
					 p3 = self.points[simplex[4]]
					 if !GoodTetOrientation(p0, p1, p2, p3)
						  simplex = [simplex[3], simplex[2], simplex[1], simplex[4]]
					 end
				end
				push!(fixed, simplex)
		  end
		  self.hulls = fixed
	 end
	 return self
end


# ////////////////////////////////////////////////////////////////////////////////
function GetBatches(self::BuildMkPol)

	sf = ToSimplicialForm(self)

	points = GLBatch(POINTS)
	lines   = GLBatch(LINES)
	triangles = GLBatch(TRIANGLES)
	
	for hull in sf.hulls
		hull_dim = length(hull)
		if hull_dim == 1
			p0 = sf.points[hull[1]];resize!(p0,3)
			push!(points.vertices, p0)
		elseif hull_dim == 2
			p0 = sf.points[hull[1]];resize!(p0,3)
			p1 = sf.points[hull[2]];resize!(p1,3)
			push!(lines.vertices.vector, p0...)
			push!(lines.vertices.vector, p1...)
		elseif hull_dim == 3
			p0 = sf.points[hull[1]];resize!(p0,3)
			p1 = sf.points[hull[2]];resize!(p1,3)
			p2 = sf.points[hull[3]];resize!(p2,3)
			n = ComputeTriangleNormal(p0, p1, p2)
			push!(triangles.vertices.vector, p0...);push!(triangles.normals.vector, n...)
			push!(triangles.vertices.vector, p1...);push!(triangles.normals.vector, n...)
			push!(triangles.vertices.vector, p2...);push!(triangles.normals.vector, n...)
		elseif hull_dim == 4
			for T in [[1, 2, 4], [1, 4, 3], [1, 3, 2], [2, 3, 4]]
					p0 = sf.points[hull[T[1]]];resize!(p0,3)
					p1 = sf.points[hull[T[2]]];resize!(p1,3)
					p2 = sf.points[hull[T[3]]];resize!(p2,3)
					n = ComputeTriangleNormal(p0, p1, p2)
					push!(triangles.vertices.vector, p0...);push!(triangles.normals.vector, n...)
					push!(triangles.vertices.vector, p1...);push!(triangles.normals.vector, n...)
					push!(triangles.vertices.vector, p2...);push!(triangles.normals.vector, n...)
			end
		else
			throw("Cannot handle geometry with dim > 3")
		end
	end
	
	ret = []
	if !isempty(points.vertices.vector)
		push!(ret, points)
	end
	if !isempty(lines.vertices.vector)
		push!(ret, lines)
	end
	if !isempty(triangles.vertices.vector)
		push!(ret, triangles)
	end
	return ret
end


# /////////////////////////////////////////////////////////////
mutable struct Hpc
	 T::MatrixNd
	 childs::Union{Vector{Hpc}, Vector{BuildMkPol}}
	 properties::Dict{Any, Any}
	 
	 # constructor
	 function Hpc(T::MatrixNd=MatrixNd(0), childs:: Union{Vector{Hpc}, Vector{BuildMkPol}}=[], properties=Dict())

		  self = new()
		  self.childs = childs
		  self.properties = properties
		  if length(childs) > 0
				Tdim = maximum([dim(child) for child in childs]) + 1
				self.T = embed(T, Tdim)
		  else
				self.T = T
		  end
		  return self
	 end
end
	 

function Base.show(io::IO, self::Hpc)
	print(io, "Hpc(", self.T, ", ", self.childs, ", ", self.properties, ")")
end

function dim(self::Hpc)
	return dim(self.T) - 1
end

function toList(self::Hpc)
	ret = []
	Tdim = dim(self) + 1
	stack = [[MatrixNd(Tdim), Dict(), self]]
	while !isempty(stack)
		T, properties, node = pop!(stack)
		if isa(node, Hpc)
				T = T * embed(node.T, Tdim)
				if !isempty(node.properties)
					properties = copy(properties)
					for (key, value) in node.properties
						properties[key] = value
					end
				end
				for child in node.childs
					push!(stack, [T, properties, child])
				end
		else
				push!(ret, [T, properties, node])
		end
	end
	return ret
end

function box(self::Hpc)
	box = BoxNd(dim(self))
	for (T, properties, obj) in toList(self)
		addPoints(box, [transformPoint(T,p) for p in obj.points])
	end
	return box
end




# //////////////////////////////////////////////////////////////////////////////////////////
function MkPol(points::Vector{Vector{Float64}}, hulls::Vector{Vector{Int}}=Vector{Vector{Int}}())
	obj=BuildMkPol(points, hulls)
	return Hpc(MatrixNd(), [obj])
end



function MkPol0()

	points=Vector{Vector{Float64}}()
	push!(points,[])
	hulls=[[1]]
	obj=BuildMkPol(points,hulls)
	return Hpc(MatrixNd(), [obj])
end

function Struct(pols::Vector{Hpc})
	return Hpc(MatrixNd(), pols)
end


# //////////////////////////////////////////////////////////////////////////////////////////
function Cube(dim::Int, From::Float64=0.0, To::Float64=1.0)

	if dim==0
		return MkPol0()
	end

	@assert dim>=1
	points = [[From],[To]]
	for I in 2:dim
		a=[ [p;From] for p in points]
		b=[ [p;To  ] for p in points]
		points =  [a;b]
	end
	return MkPol(points)
end

# //////////////////////////////////////////////////////////////////////////////////////////
function Simplex(dim::Int)

	if dim==0
		return MkPol0()
	end

	points = [ [0.0 for _ in 1:dim] ]
	for I in 1:dim
		point=[0.0 for _ in 1:dim]
		point[I]=1.0
		push!(points,point)
	end
	return MkPol(points)
end

# //////////////////////////////////////////////////////////////////////////////////////////
function Join(pols::Vector{Hpc})
	points = Vector{Vector{Float64}}()
	for (T, properties, obj) in toList(Hpc(MatrixNd(), pols))
		append!(points, [transformPoint(T,p) for p in obj.points])
	end
	return MkPol(points)
end

# //////////////////////////////////////////////////////////////////////////////////////////
function Quote(sequence::Vector{Float64})
	pos = 0.0
	points = [[pos]]
	hulls = Vector{Vector{Int}}()
	for value in sequence
		next = pos + abs(value)
		push!(points, [next])
		if value >= 0
				push!(hulls, [length(points)-1, length(points)])
		end
		pos = next
	end
	return MkPol(points, hulls)
end


# //////////////////////////////////////////////////////////////////////////////////////////
function Transform(self::Hpc, T::MatrixNd)
	return Hpc(T, [self])
end

function Translate(self::Hpc, vt::Vector{Float64})
	return Hpc(translate(vt), [self])
end

function Scale(self::Hpc, vs::Vector{Float64})
	return Hpc(scale(vs), [self])
end

function Rotate(self::Hpc, i::Int, j::Int, angle::Float64)
	return Hpc(rotate(i, j, angle), [self])
end

# //////////////////////////////////////////////////////////////////////////////////////////
function Power(a::Hpc, b::Hpc)
	childs = Vector{Hpc}()
	for (T2, properties2, obj2) in toList(b)
		for (T1, properties1, obj1) in toList(a)
				points = Vector{Vector{Float64}}()
				for py in obj2.points
					for px in obj1.points
						push!(points, [px;py])
					end
				end
				hulls = Vector{Vector{Int}}()
				nx, ny = length(obj1.points), length(obj2.points)
				for hy in obj2.hulls
					for hx in obj1.hulls
						hull=Vector{Int}()
						for A in hx
							for B in hy
								push!(hull,1+((B-1)*nx + (A-1)))
							end
						end
						push!(hulls,hull)
					end
				end
				T = adjoin(T1, T2)
				push!(childs, Hpc(T, [BuildMkPol(points, hulls)]))
		end
	end
	return Hpc(MatrixNd(), childs)
end

# //////////////////////////////////////////////////////////////////////////////////////////
function UkPol(self::Hpc)
	points = Vector{Vector{Float64}}()
	hulls  = Vector{Vector{Int}}()
	for (T, properties, obj) in toList(self)
		O = length(points)

		for p in obj.points
			push!(points,transformPoint(T,p))
		end

		for hull in obj.hulls
			new_hull=[O + idx for idx in hull]
			push!(hulls,new_hull)
		end

	end

	# incredible that Julia converts to a vector of vector of float automatically
	# return [points, hulls]

	ret=Vector{Any}()
	push!(ret,points)
	push!(ret,hulls)
	return ret


end

# //////////////////////////////////////////////////////////////////////////////////////////
function GetBatches(self::Hpc)
	batches = Vector{GLBatch}()
	for (T, properties, obj) in toList(self)
		T = embed(T, 4)
		T4d=Matrix4d(
			T[2,2], T[2,3], T[2,4],   T[2,1],  # homo should be last
			T[3,2], T[3,3], T[3,4],   T[3,1],
			T[4,2], T[4,3], T[4,4],   T[4,1],

			T[1,2], T[1,3], T[1,4],   T[1,1]
		)
		for batch in GetBatches(obj)
				prependTransformation(batch, T4d)
				# writeProperties(batch, properties) TODO: loosing properties
				push!(batches, batch)
		end
	end
	return batches
end


# //////////////////////////////////////////////////////////////////////////////////////////
function View(self::Hpc, title::String="Plasm.jl")
	batches=[GLAxis(Point3d(0,0,0),Point3d(+2,2,2));GetBatches(self)]
	GLView(batches, title)
end

# //////////////////////////////////////////////////////////////////////////////////////////
function MapFn(self::Hpc, fn)
	childs = Vector{Hpc}()
	for (T, properties, obj) in toList(self)
		sf = ToSimplicialForm(obj)
		points = [fn(transformPoint(T,p)) for p in sf.points]
		hulls = sf.hulls
		push!(childs, Hpc(MatrixNd(), [BuildMkPol(points, hulls)], properties))
	end
	ret = Hpc(MatrixNd(), childs)
	return ret
end

# //////////////////////////////////////////////////////////////////////////////////////////
function ToBoundaryForm(self::Hpc)
	DB = Dict{Vector{Float64},Int64}()
	POINTS=Vector{Vector{Float64}}()
	FACES = Vector{Vector{Int}}()

	for (T, properties, obj) in toList(self)
		sf = ToSimplicialForm(obj)
		points, hulls = [transformPoint(T,p) for p in sf.points], sf.hulls
		pdim = length(points[1])
		mapped = Dict{Int64,Int64}()


		for P in 1:length(points)
				point = points[P]
				@assert point isa Vector{Float64}
				idx=get(DB,point,0)
				if idx==0
					push!(POINTS,point)
					DB[point]=length(POINTS)
				end
				mapped[P]=DB[point]
		end

		for hull in hulls
				bfaces = []
				if length(hull) < (pdim + 1)
					bfaces = [collect(1:length(hull))]
				else
					@assert  length(hull) == (pdim + 1) # it shouls be a simplex
					if pdim == 0
						bfaces = [[1]]
					elseif pdim == 1
						bfaces = [[1], [2]]
					elseif pdim == 2
						bfaces = [[1, 2], [2, 3], [3, 1]]
					elseif pdim == 3
						bfaces = [[1, 2, 4], [1, 4, 3], [1, 3, 2], [2, 3, 4]]
					else
						error("not supported")
					end
				end
				for face in bfaces
					push!(FACES, [mapped[hull[it]] for it in face])
				end
		end
	end

	num_occurrence = Dict{Vector{Int64},Int}()
	for face in FACES
		Key = sort(face)
		num_occurrence[Key] = get(num_occurrence,Key,0) + 1
	end
 
	FACES=[face for face in FACES if num_occurrence[sort(face)] == 1]
	return MkPol(POINTS,FACES)
end


# //////////////////////////////////////////////////////////////////////////////////
# https://github.com/scipy/scipy/blob/main/scipy/spatial/_qhull.pyx
# ConvexHull class triangulate the resutls (see required_optioins) 
# so I am creating a new class to override the problem


# //////////////////////////////////////////////////////////////////////////////////////////
function ToLARForm(self::Hpc)

	py"""

	import numpy as np
	import math
	from scipy.spatial._qhull import _Qhull,_QhullUser
	
	# //////////////////////////////////////////////////////////////////////////////////
	class MyConvexHull(_QhullUser):
	
		# constructor
		def __init__(self, user_points, verbose=False):
			user_points = np.ascontiguousarray(user_points, dtype=np.double)
			qhull_options = "i Qs Pp".encode('latin1')
			qhull = _Qhull(b"i", user_points, qhull_options, required_options=None, incremental=True) # this is the line I need to change
			_QhullUser.__init__(self, qhull, incremental=True)
	
			self.pdim=len(user_points[0])
	
			self.points=self._points
			self.facets,self.normals=qhull.get_hull_facets()
			self.close()
	
			# in 2d the normal is pointing up (i.e Z>=0)
			if self.pdim==2:
				self.normals=[[0.0,0.0,1.0]]
		
			self.removeUnusedPoints()
			self.findEdges()
			self.findFacetLoops()
			self.correctOrientation()
	
			if verbose:
				print("Points", len(self.points))
				for P, point in enumerate(self.points): 
					print(P,point)
	
				print("Facets", len(self.facets))
				for F,facet in enumerate(self.facets): 
					print(facet, self.normals[F])
	
		# GetPlane
		@staticmethod
		def GetPlane(p1,p2,p3): 
			x1,y1,z1=p1
			x2,y2,z2=p2
			x3,y3,z3=p3
			a1 = x2 - x1
			b1 = y2 - y1
			c1 = z2 - z1
			a2 = x3 - x1
			b2 = y3 - y1
			c2 = z3 - z1
			a = b1 * c2 - b2 * c1
			b = a2 * c1 - a1 * c2
			c = a1 * b2 - b1 * a2
			d = (- a * x1 - b * y1 - c * z1)
			m=math.sqrt(a*a+b*b+c*c)
			if not m: m=1.0
			return [a/m,b/m,c/m,d/m]
	
		#  removeUnusedPoints
		def removeUnusedPoints(self):
			reindex,new_points={},[]
			for F,facet in enumerate(self.facets):
				for idx in facet:
					if not idx in reindex:
						new_index=len(reindex)
						reindex[idx]=new_index
						new_points.append(self.points[idx])
				self.facets[F]=[reindex[idx] for idx in facet]
			self.points=new_points
	
		# findEdges
		def findEdges(self):
			if self.pdim==2:
	
				# only one face
				self.edges={0: self.facets} 
	
			elif self.pdim==3:
	
				# find edges as the intersection of faces
				if True:
					nfaces=len(self.facets)
					self.edges={I:[] for I in range(nfaces)}
					for A in range(nfaces):
						for B in range(A+1,nfaces):
							edge=list(set.intersection(set(self.facets[A]),set(self.facets[B])))
							if edge:
								self.edges[A].append(edge)
								self.edges[B].append(edge)
	
		# findFacetLoops
		def findFacetLoops(self):
			new_facets=[]
			for F in self.edges:
				faces=self.facets[F]
				todo=self.edges[F]
				loop=[]
				(A,B),todo=todo[0],todo[1:]
				loop.append(A)
				while todo:
					found=False
					for I,(a,b) in enumerate(todo):
						if a==B or b==B:
							loop.append(B)
							B=b if a==B else a
							found=True
							del todo[I]
							break
					assert(found)
				new_facets.append(loop)
			self.facets=new_facets
	
		# DotProduct
		@staticmethod
		def DotProduct(a,b):
			return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
	
		# correctOrientation
		def correctOrientation(self):
			for F,face in enumerate(self.facets):
				p1,p2,p3=[self.points[idx] for idx in self.facets[F][0:3]]
				if self.pdim==2:
					p1,p2,p3=list(p1)+[0.0],list(p2)+[0.0],list(p3)+[0.0]
				plane=MyConvexHull.GetPlane(p1,p2,p3)
				if MyConvexHull.DotProduct(plane,self.normals[F])<0:
					self.facets[F].reverse()
	
	def GetLARConvexHull(user_points, verbose=False):
		lar=MyConvexHull(user_points,verbose=False)
		return [
			lar.points, 
			[[int(it+1) for it in facet] for facet in lar.facets] # 0->1 index
		]
	
	"""

	ret=BuildMkPol()
	for (T, properties, obj) in toList(self)

		@assert obj isa BuildMkPol

		# automatic filter useless obj
		if isempty(obj.points) || isempty(obj.hulls) || dim(obj) <= 1
			continue
		end

		# for each hull create a LAR hull (i.e. polygonal faces)
		for hull in obj.hulls
			points=[obj.points[idx] for idx in hull]
			points, facets = py"GetLARConvexHull"(points)

			#println(typeof(points),points)
			#println(typeof(facets),facets)

			if facets isa Matrix{Int64}
				nrows,ncols=size(facets)
				converted=Vector{Vector{Int64}}()
				for R in 1:nrows
					push!(converted, facets[R,:])
				end
				facets=converted
			end

			@assert points isa Vector{Vector{Float64}}
			@assert facets isa Vector{Vector{Int64}}
			
			# add the transformed points
			points = [transformPoint(T,p) for p in points]
			mapped = addPoints(ret, points)

			# add the hull
			push!(ret.hulls, [mapped[it] for it in 1:length(points)])

			# add the facets
			for facet in facets
				push!(ret.facets, [mapped[it] for it in facet])
			end
		end
	end
	return Hpc(MatrixNd(), [ret])
end

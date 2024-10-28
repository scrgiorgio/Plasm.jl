const BSP_EPSILON=1e-5

# /////////////////////////////////////////////////////////////////////////////
struct BspFace

	vertices::Vector{PointsNd}
	plane::PointNd

  # constructor
	function BspFace(vertices, plane=nothing)

		# automatic computation of plane
		if isnothing(plane)
			dim = length(vertices[1])
			if dim == 2
				@assert length(vertices) == 2
				x1, y1 = vertices[1][1], vertices[1][2]
				x2, y2 = vertices[2][1], vertices[2][2]
				normal = [(y2 - y1), -(x2 - x1)]
			elseif dim == 3
				normal = ComputeTriangleNormal(vertices[1], vertices[2], vertices[3])
			else
				error("todo")
			end

			w = sum([normal[i] * vertices[1][i] for i in 1:dim])
			plane = [normal; -w]
		end
		new(vertices, plane)
	end

end

# //////////////////////////////////////////////////////////////
function split_edge(plane::PointNd, vi::PointNd, vj::PointNd)
	dim = length(vi)
	valuev1 = abs(plane[end] + sum([plane[i] * vi[i] for i in 1:dim]))
	valuev2 = abs(plane[end] + sum([plane[i] * vj[i] for i in 1:dim]))
	alpha = 1.0 / (valuev1 + valuev2)
	beta = valuev1 * alpha
	alpha = alpha * valuev2
	return [alpha * vi[i] + beta * vj[i] for i in 1:dim]
end

# //////////////////////////////////////////////////////////////
function split_face(face::BspFace, plane::PointNd)
	dim = length(plane) - 1
	COPLANAR, ABOVE, BELOW, SPANNING = 0, 1, 2, 3
	ptype, types = COPLANAR, []
	for v in face.vertices
		@assert length(v) == dim
		t = plane[end] + sum([plane[i] * v[i] for i in 1:dim])
		if t < -BSP_EPSILON
			type = BELOW
		elseif t > BSP_EPSILON
			type = ABOVE
		else
			type = COPLANAR
		end
		ptype |= type
		push!(types, type)
	end

  # return [<epsilon,coplanar below, coplanar above,>epsilon]
	if ptype == BELOW
		return [face, nothing, nothing, nothing]
	end
	
  if ptype == ABOVE
		return [nothing, nothing, nothing, face]
	end
	
  if ptype == COPLANAR
		if sum([plane[i] * face.plane[i] for i in 1:dim+1]) > 0
			return [nothing, nothing, face, nothing]
		else
			return [nothing, face, nothing, nothing]
		end
	end

  # need to cut
	@assert ptype == SPANNING
	b, f = [], []

  # split and edge
	if dim == 2
		@assert length(face.vertices) == 2
		ti, tj = types[1], types[2]
		vi, vj = face.vertices[1], face.vertices[2]
		if ti != BELOW
			push!(f, vi)
		end
		if ti != ABOVE
			push!(b, vi)
		end
		if tj != BELOW
			push!(f, vj)
		end
		if tj != ABOVE
			push!(b, vj)
		end
		if (ti | tj) == SPANNING
			v = split_edge(plane, vi, vj)
			push!(b, v)
			push!(f, v)
		end

  # split a face 
	elseif dim == 3
		for i in 1:length(face.vertices)
			j = (i + 1) % length(face.vertices)
			ti, tj = types[i], types[j]
			vi, vj = face.vertices[i], face.vertices[j]
			if ti != BELOW
				push!(f, vi)
			end
			if ti != ABOVE
				push!(b, vi)
			end
			if (ti | tj) == SPANNING
				v = split_edge(plane, vi, vj)
				push!(b, v)
				push!(f, v)
			end
		end
	else
		error("not supported")
	end
	@assert length(b) >= dim && length(f) >= dim
	return [
    BspFace(b, face.plane), 
    nothing, 
    nothing, 
    BspFace(f, face.plane)
  ]
end

# //////////////////////////////////////////////////////////////////////////////////////
struct Bsp

	plane::PointNd
	faces::Vector{BspFace}
	below::Bsp
	above::Bsp

  # constructor
	function Bsp()
		new()
	end

end

# ////////////////////////////////////////////////////////////////////////
function bsp_get_faces(self::Bsp)::Vector{BspFace}
	if isnothing(self) return [] end
	return [self.faces; bsp_get_faces(self.below); bsp_get_faces(self.above)]
end

# ////////////////////////////////////////////////////////////////////////
function bsp_insert_faces(self::Bsp, faces::Vector{BspFace})
	if length(faces) == 0
		return self
	end

  # this will be the root
	if isnothing(self.plane)
		@assert isnothing(self.below)
    @assert isnothing(self.above)
		self.plane = faces[1].plane
	end

	below, above = [], []
	for face in faces
		b, cb, ca, a = bsp_split_face(face, self.plane)
		if !isnothing(b ) push!(below,       b) end
		if !isnothing(cb) push!(self.faces, cb) end
		if !isnothing(ca) push!(self.faces, ca) end
		if !isnothing(a ) push!(above,       a) end
	end

	if !isempty(above)
		if isnothing(self.above) self.above = Bsp() end
		bsp_insert_faces(self.above, above)
	end

	if !isempty(below)
		if isnthing(self.below) self.below = Bsp() end
		bsp_insert_faces(self.below, below)
	end
end

# /////////////////////////////////////////////////////////////////
function bsp_fragment_faces(self::Bsp, faces::Vector{BspFace})::Vector{BspFace}
	if isnothing(self.plane) 
    return faces 
  end

	below, above = [], []
	for face in faces
		b, cb, ca, a = split_face(face, self.plane)
		if !isnothing( b) push!(below,  b) end
		if !isnothing(cb) push!(below, cb) end
		if !isnothing(ca) push!(above, ca) end
		if !isnothing( a) push!(above,  a) end
	end
	below = !isnothing(self.below) ? bsp_fragment_faces(self.below, below) : []
	above = !isnothing(self.above) ? bsp_fragment_faces(self.above, above) : []
	return above + below
end

# /////////////////////////////////////////////////////////////////
function bsp_clip_to(a::Bsp, b::Bsp)
	a.faces = bsp_fragment_faces(b, a.faces)
	if !isnothing(a.below) bsp_clip_to(a.below, b) end
	if !isnothing(a.above) bsp_clip_to(a.above, b) end
end

# /////////////////////////////////////////////////////////////////
function bsp_complement(bsp::Bsp)::Bsp
	if isnothing(bsp) return nothing end
	ret = Bsp()

  # reverse plane
	if !isnothing(bsp.plane)
		ret.plane = [-1 * it for it in bsp.plane]
	end

  # reverse faces
	for p in bsp.faces
		new_p = BspFace(reverse(p.vertices), [-1 * c for c in p.plane])
		push!(ret.faces, new_p)
	end
	ret.below = BspComplement(bsp.above)
	ret.above = BspComplement(bsp.below)
	return ret
end

# /////////////////////////////////////////////////////////////////////
function bsp_union(a::Bsp, b::Bsp)::Bsp
	bsp_clip_to(a, b)
	bsp_clip_to(b, a)
	b = BspComplement(b)
	bsp_clip_to(b, a)
	b = bsp_complement(b)
	bsp_insert_faces(a, bsp_get_faces(b))
	return a
end

# /////////////////////////////////////////////////////////////////////
function bsp_intersection(a::Bsp, b::Bsp)::Bsp
	return bsp_complement(
    bsp_union(
      bsp_complement(a), 
      bsp_complement(b)
    )
  )
end

# /////////////////////////////////////////////////////////////////////
function bsp_difference(a::Bsp, b::Bsp)
	return bsp_complement(
    bsp_union(
      bsp_complement(a), 
      b))
end

# /////////////////////////////////////////////////////////////////////
function bsp_xor(a::Bsp, b::Bsp)
	return bsp_union(
    bsp_intersection(
      a, 
      bsp_complement(b)
    ), 
    bsp_intersection(
      bsp_complement(a), 
      b
    )
  )
end


# /////////////////////////////////////////////////////////////////////
function bsp_from_hpc(hpc::Hpc)::Bsp

	geo = ToGeometry(hpc)

  points=geo.points
  faces=geo.faces

  ret = Bsp()
	bsp_insert_faces(ret, faces)
	return ret
end

# /////////////////////////////////////////////////////////////////////
function bsp_to_hpc(self::Bsp)::Hpc
  batches = []
	faces = bsp_get_faces(self)
	dim = self.plane !== nothing ? length(self.plane) - 1 : 0

	if dim == 0
		return Hpc()
	end

  points, hulls = [], []

  if dim == 1
    for face in faces
			@assert length(face.vertices) == 1
			N = length(points)
			points += face.vertices
			push!(hulls, collect(N:N+length(face.vertices)-1))      
    end
    return MkPol(points, hulls)
  end

  if dim == 2
    for face in faces
			@assert length(face.vertices) == 2
			N = length(points)
			points += face.vertices
			push!(hulls, collect(N:N+length(face.vertices)-1))
    end
    return MkPol(points, hulls)
  end

  if dim == 3
    for face in faces
			@assert length(face.vertices) >= 3
			for I in 2:length(face.vertices)-1
				N = length(points)
				points += [face.vertices[1], face.vertices[I], face.vertices[I+1]]
				push!(hulls, collect(N:N+2))
			end
      return MkPol(points, hulls)

    end
  end
  error("not supported")

end

# /////////////////////////////////////////////////////////////////////
function bsp_union(objs::Vector{Hpc})::Hpc
	objs = [bsp_from_hpc(obj) for obj in objs]
	res = nothing
	for it in objs
		res = isnothing(res) ? it : bsp_union(res, it)
	end
	return bsp_to_hpc(res)
end

# /////////////////////////////////////////////////////////////////////
function bsp_intersection(objs::Vector{Hpc})::Hpc
	objs = [bsp_from_hpc(obj) for obj in objs]
	res = nothing
	for it in objs
		res = isnothing(res) ? it : bsp_intersection(res, it)
	end
	return bsp_to_hpc(res)
end

# /////////////////////////////////////////////////////////////////////
function bsp_difference(objs::Vector{Hpc})::Hpc
	objs = [bsp_from_hpc(obj) for obj in objs]
	res = nothing
	for it in objs
		res = isnothing(res) ? it : bsp_difference(res, objs[I])
	end
	return bsp_to_hpc(res)
end

# /////////////////////////////////////////////////////////////////////
function bsp_xor(objs::Vector{Hpc})::Hpc
	objs = [bsp_from_hpc(obj) for obj in objs]
	res = nothing
	for it in objs
		res = isnothing(res) ? it : bsp_xor(res, objs[I])
	end
	return bsp_to_hpc(res)
end



orthogonal_axis=[[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0] ]
export orthogonal_axis

# ////////////////////////////////////////////////////////////////////////
""" needed for normals """
function normalized(value::Vector{Float64})::Vector{Float64}
  return value/LinearAlgebra.norm(value)
end
export normalized

# ////////////////////////////////////////////////////////////////////////
""" note: plane normal should be normalized """
function plane_point_distance(plane::Vector{Float64},point::Vector{Float64})
  A,B,C,D=plane
  x,y,z=point
  return A*x+B*y+C*z+D
end
export random_plane

# ////////////////////////////////////////////////////////////////////////
""" generate random plane"""
function random_plane()
  origin = rand(3)
  normal = normalized(rand(3))
  A,B,C=normal
  D=-(A*origin[1]+B*origin[2]+C*origin[3])
  plane=[A,B,C,D]
  @assert(abs(plane_point_distance(plane,origin))<Plasm.LAR_DEFAULT_ERR)
  return plane
end
export random_plane


# ////////////////////////////////////////////////////////////////////////
""" generate random points on plane """
function random_points_on_plane(plane::Vector{Float64}; num_vertices::Int)
  normal=[plane[1],plane[2],plane[3]]
  points=[]
  for I in 1:num_vertices
    point=rand(3)
    dist=plane_point_distance(plane,point)
    point=point - dist*normal
    @assert(abs(plane_point_distance(plane,point))<Plasm.LAR_DEFAULT_ERR)
    push!(points, point)
  end
  return stack(points,dims=2) 

end
export random_points_on_plane

# ////////////////////////////////////////////////////////////////////////
function face_coordinate_system(V::Points)
  center = compute_centroid(V)
  V=stack([point-center for point in eachcol(V)],dims=2)
  __,__,normal=[normalize(it) for it in eachrow(svd(BYROW(V)).Vt)]
  v=normalized(cross(normal,orthogonal_axis[argmin(normal)]))
  u=normalized(cross(v,normal))
  return center,u,v,normal
end
export face_coordinate_system

orthogonal_axis=[[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0] ]
export orthogonal_axis

# ////////////////////////////////////////////////////////////////////////
""" needed for normals """
function normalized(value::Vector{Float64})::Vector{Float64}
  return value/LinearAlgebra.norm(value)
end
export normalized

# ////////////////////////////////////////////////////////////////////////
""" create plane """
function plane_create(normal::Vector{Float64},point::Vector{Float64})::Vector{Float64}
  A,B,C=normalized(normal)
  D=-dot(normal,point)
  return [A,B,C,D]
end
export plane_create


# ////////////////////////////////////////////////////////////////////////
function plane_create(V::Points) ::Vector{Float64}
  center = compute_centroid(V) 
  V=stack([point-center for point in eachcol(V)],dims=2)
  __,__,normal=[normalize(it) for it in eachrow(svd(BYROW(V)).Vt)]
  return plane_create(normal, center)
end

# ////////////////////////////////////////////////////////////////////////
function plane_get_normal(plane::Vector{Float64})
  A,B,C,D=plane
  return [A,B,C]
end
export plane_get_normal

# ////////////////////////////////////////////////////////////////////////
""" note: plane normal should be normalized """
function plane_point_distance(plane::Vector{Float64},point::Vector{Float64})
  A,B,C,D=plane
  x,y,z=point
  return A*x+B*y+C*z+D
end
export plane_point_distance

# ////////////////////////////////////////////////////////////
""" get ray and plane intersection """
function plane_ray_intersection(ray_origin::Vector{Float64}, ray_dir::Vector{Float64}, plane::Vector{Float64})

  ray_dir=normalized(ray_dir)

  # see https://education.siggraph.org/static/HyperGraph/raytrace/rayplane_intersection.htm
  denom=dot(plane_get_normal(plane),ray_dir)

  # just because t would go to infinite...
  if abs(denom)<LAR_DEFAULT_ERR
    return nothing
  else
    t= -plane_point_distance(plane,ray_origin)/denom
    if t<=0
      return nothing # in below space, ignore
    else
      return ray_origin+t*ray_dir
    end
  end

end
export plane_ray_intersection


# ////////////////////////////////////////////////////////////////////////
""" generate random points on plane """
function plane_random_points(plane::Vector{Float64}; num_vertices::Int)
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
export plane_random_points



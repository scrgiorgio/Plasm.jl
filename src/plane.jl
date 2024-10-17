
orthogonal_axis=[[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0] ]
export orthogonal_axis


# ////////////////////////////////////////////////////////////////////////
""" needed for normals """
function normalized(value::PointNd)::PointNd
  return value/LinearAlgebra.norm(value)
end
export normalized

# ////////////////////////////////////////////////////////////////////////
""" create plane """
function plane_create(normal::PointNd, point::PointNd)::PointNd
  A,B,C=normalized(normal)
  D=-dot(normal,point)
  return [A,B,C,D]
end
export plane_create


# ////////////////////////////////////////////////////////////////////////
function plane_create(V::Points) ::PointNd
  center = compute_centroid(V) 
  V=stack([point-center for point in eachcol(V)],dims=2)
  __,__,normal=[normalize(it) for it in eachrow(svd(BYROW(V)).Vt)]
  return plane_create(normal, center)
end

# ////////////////////////////////////////////////////////////////////////
function plane_get_normal(plane::PointNd)
  A,B,C,D=plane
  return [A,B,C]
end
export plane_get_normal

# ////////////////////////////////////////////////////////////////////////
""" note: plane normal should be normalized """
function plane_point_distance(plane::PointNd,point::PointNd)
  A,B,C,D=plane
  x,y,z=point
  return A*x+B*y+C*z+D
end
export plane_point_distance

# ////////////////////////////////////////////////////////////
""" get ray and plane intersection """
function plane_ray_intersection(ray_origin::PointNd, ray_dir::PointNd, plane::PointNd)

  ray_dir=normalized(ray_dir)

  # see https://education.siggraph.org/static/HyperGraph/raytrace/rayplane_intersection.htm
  denom=dot(plane_get_normal(plane),ray_dir)

  # just because t would go to infinite...
  if abs(denom)<LAR_DEFAULT_ERR
    return nothing,nothing
  else
    t= -plane_point_distance(plane,ray_origin)/denom
    if t<=0
      return nothing, nothing# in below space, ignore
    else
      return ray_origin+t*ray_dir,t
    end
  end

end
export plane_ray_intersection


# ////////////////////////////////////////////////////////////////////////
""" generate random points on plane """
function plane_random_points(plane::PointNd; num_vertices::Int)
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


# ///////////////////////////////////////////////////////////////
function project_points3d(V::Points; double_check=false)

  C = compute_centroid(V)

  # points reference frame
  Z=plane_get_normal(plane_create(V))
  Y=normalized(cross(Z,orthogonal_axis[argmin([abs(it) for it in Z])]))
  X=normalized(cross(Y,Z))

  # if I go from 1,0,0  0,1,0  0,0,1 to XYZ I need a Tdir, Tinv is the transposed
  Tdir=[
    X[1] Y[1] Z[1] ; 
    X[2] Y[2] Z[2] ; 
    X[3] Y[3] Z[3]
  ]
  
  Tinv=transpose(Tdir)

  function internal_projector(V::Points; inverse=false, keep_z=false)

    if inverse
      ret=[[x,y,0.0] for (x,y) in eachcol(V)]
      ret=[((Tdir * p) + C) for p in ret]
      return hcat(ret...) # this returns by-col points2d
    else

      points3d=[p for p in eachcol(V)]
      ret=[(Tinv * (p-C)) for p in points3d]

      # try if applying the direct will return the same points
      if double_check
        checking=[(Tdir * p)+C for p in ret]
        for (a,b) in zip(points3d, checking)
          @assert(vertex_fuzzy_equals(a,b))
        end
      end

      # remove Z coordinate
      if !keep_z
        ret=[p[1:2] for p in ret]
      end

      # this returns by-col
      return hcat(ret...)
    end

  end

  return internal_projector

end
export project_points3d

using Test
using Plasm, LinearAlgebra

@testset "plane.jl" begin

  @testset "face coordinate system" begin
    for I in 1:10
      plane=plane_create(rand(3),rand(3))
      V=plane_random_points(plane,num_vertices=100)
      center,v1,v2,v3=face_coordinate_system(V)
      @assert(vertex_fuzzy_equals(v3,plane[1:3]) || vertex_fuzzy_equals(-v3,plane[1:3]))
    end
  end



  @testset "ray" begin
  
    ray_origin= [0.0,0.0,-1.0]
    ray_dir=normalized([0.0, 0.0, 1.0])

    plane_normal=normalized([0.0, 0.0, 1.0])
    plane_normal_inv=-1*plane_normal
    plane_point=[0.0,0.0,2.0]

    # it does not matter the orientation of the plane
    point=plane_ray_intersection(ray_origin, +1*ray_dir , plane_create(plane_normal,plane_point))
    @assert(point==[0.0,0.0,2.0])

    point=plane_ray_intersection(ray_origin, +1*ray_dir, plane_create(plane_normal_inv,plane_point))
    @assert(point==[0.0,0.0,2.0])

    # if rays is going in the opposite direction, techincally we don't have any intersection
    point=plane_ray_intersection(ray_origin, -1*ray_dir, plane_create(plane_normal,plane_point))
    @assert(isnothing(point))

    point=plane_ray_intersection(ray_origin, -1*ray_dir, plane_create(plane_normal_inv,plane_point))
    @assert(isnothing(point))

  end

end;

using Test
using Plasm, LinearAlgebra

@testset "plane.jl" begin

  @testset "face coordinate system" begin
    for I in 1:10
      plane=plane_create(rand(3),rand(3))
      V=plane_random_points(plane,num_vertices=100)
      normal      =plane_get_normal(plane)
      normal_check=plane_get_normal(plane_create(V))
      @test(vertex_fuzzy_equals(normal,plane[1:3]) || vertex_fuzzy_equals(-normal,plane[1:3]))
    end
  end

  @testset "ray" begin
  
    ray_origin= [0.0,0.0,-1.0]
    ray_dir=normalized([0.0, 0.0, 1.0])

    plane_normal=normalized([0.0, 0.0, 1.0])
    plane_normal_inv=-1*plane_normal
    plane_point=[0.0,0.0,2.0]

    # it does not matter the orientation of the plane
    point, t=plane_ray_intersection(ray_origin, +1*ray_dir , plane_create(plane_normal,plane_point))
    @test(point==[0.0,0.0,2.0])

    point, t=plane_ray_intersection(ray_origin, +1*ray_dir, plane_create(plane_normal_inv,plane_point))
    @test(point==[0.0,0.0,2.0])

    # if rays is going in the opposite direction, techincally we don't have any intersection
    point, t=plane_ray_intersection(ray_origin, -1*ray_dir, plane_create(plane_normal,plane_point))
    @test(isnothing(point))

    point, t=plane_ray_intersection(ray_origin, -1*ray_dir, plane_create(plane_normal_inv,plane_point))
    @test(isnothing(point))

  end

  @testset "projector" begin

    x,y,z=1.0, 2.0, 3.0
  
    V=BYCOL([
      x-1.0  y  z-1.0;
      x+1.0  y  z-1.0;
      x+1.0  y  z+1.0;
      x-1.0  y  z+1.0
    ])
    projector=project_points3d(V;double_check=true)
    @assert(size(projector(V),2)==4)
  
    V=BYCOL([
      x-1.0  y-1.0  z;
      x+1.0  y-1.0  z;
      x+1.0  y+1.0  z;
      x-1.0  y+1.0  z
    ])
    projector=project_points3d(V;double_check=true)
    @assert(size(projector(V),2)==4)
  
    V=BYCOL([
      x  y-1.0  z-1.0;
      x  y+1.0  z-1.0;
      x  y+1.0  z+1.0;
      x  y-1.0  z+1.0
    ])
    projector=project_points3d(V;double_check=true)
    @assert(size(projector(V),2)==4)

    for I in 1:10
      plane=plane_create(rand(3),rand(3))
      points=plane_random_points(plane, num_vertices=100)
      local projector=project_points3d(points; double_check=true)
      projector(points)
    end
  
  end


end;
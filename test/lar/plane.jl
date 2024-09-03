
using Test
using Plasm, LinearAlgebra

@testset "plane.jl" begin

  @testset "face coordinate system" begin
    for I in 1:10
      plane=random_plane()
      V=random_points_on_plane(plane,num_vertices=100)
      center,v1,v2,v3=face_coordinate_system(V)
      @assert(vertex_fuzzy_equals(v3,plane[1:3]) || vertex_fuzzy_equals(-v3,plane[1:3]))
    end
  end

end;
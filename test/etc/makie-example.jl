using GLMakie
using GeometryBasics

scene = Scene()

cam = Makie.Camera3D(scene, projectiontype = Makie.Orthographic, clipping_mode= :static)

update_cam!(scene, [3.0,3.0,3.0], [0.0,0.0,0.0],  [0.0, 0.0, 1.0])



# @show(cam)


# lines
if false
  linesegments([
      Point2f0(0, 0) => Point2f0(1, 0);
      Point2f0(1, 0) => Point2f0(1, 1);
      Point2f0(1, 1) => Point2f0(0, 1);
      Point2f0(0, 1) => Point2f0(1, 0);
    ], 
    color = :red, 
    linewidth = 2)
end

# text
text!(scene, "HELLO", position = Vec3f0(0.5,0.5,0.5), space = :data)

# triangles
a=mesh!(scene,
  [
    0.0 0.0 0.0;
    1.0 0.0 0.0;
    1.0 1.0 0.0;
    0.0 1.0 0.0
  ], 
  [1 2 3;3 4 1;], 
  # normal=[0.0 0.0 1.0;0.0 0.0 1.0;0.0 0.0 1.0;0.0 0.0 1.0;],
  color = [:red, :green, :blue, :orange])

# scene
gl_screen = display(scene)
wait(gl_screen)


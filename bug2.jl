using Plasm
using Random

Random.seed!(0)

assembly = STRUCT(
  CUBE(1), 
  T(1,2,3)(.5,.5,.5), 
  CUBE(1)
)

lar=LAR(assembly)
VIEWCOMPLEX(lar)

arrangement = ARRANGE3D(lar)
VIEWCOMPLEX(arrangement, explode=[1.4,1.4,1.4])

atoms=[]
for face_indices in arrangement.C[:CF]
  atom=SELECT(arrangement, face_indices)
  VIEWCOMPLEX(atom, explode=[1.4,1.4,1.4])
  push!(atoms,atom)
end

aaa()

internal_points=[ [internal_point,ray_dir] for (internal_point, ray_dir) in atoms] 

for (atom, (ray_origin, ray_dir)) in zip(atoms, internal_points)

  begin
    points = GLBatch(POINTS)
    points.point_size = 3
    points.point_color = RED
    append(points.vertices.vector, ray_origin)
  end

  begin
    lines = GLBatch(LINES)
    lines.line_width = 3
    lines.point_color = YELLOW
    append(lines.vertices.vector, (ray_origin+0.0*ray_dir))
    append(lines.vertices.vector, (ray_origin+1.0*ray_dir))
  end

  batches=BATCHES(atom)
  append!(batches,[points,lines])
  VIEWBATCHES(batches)
end


# bool3d(assembly, arrangement, CF)

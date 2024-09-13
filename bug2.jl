using Plasm
using Random
using LinearAlgebra

Random.seed!(0)

assembly = STRUCT(
  CUBE(1), 
  T(1,2,3)(.5,.5,.5), 
  CUBE(1)
)

lar=LAR(assembly)
# VIEWCOMPLEX(lar)

arrangement = ARRANGE3D(lar)
# VIEWCOMPLEX(arrangement, explode=[1.4,1.4,1.4])

atom_faces=arrangement.C[:CF]

atoms=[]
for sel in atom_faces
  atom=SELECT(arrangement, sel)
  push!(atoms,atom)
  # VIEWCOMPLEX(atom, explode=[1.4,1.4,1.4])
end

# remove outer atoms (i.e. the one that contains all other atoms) by using diagonal measurement
begin
  diags =[LinearAlgebra.norm(b[2] - b[1]) for b in [lar_bounding_box(atom, only_used_vertices=true) for atom in atoms]]
  __outer, outer_index = findmax(diags)
  deleteat!(atoms,      outer_index)
  deleteat!(atom_faces, outer_index)
end

internal_points=[find_internal_point(atom) for atom in atoms] 

input_args = [LAR(it) for it in TOPOS(assembly)]
bool_matrix=zeros(Bool,length(atoms),length(input_args))

for (A,(atom, (internal_point, ray_dir, distances))) in enumerate(zip(atoms, internal_points))
  @show(internal_point, ray_dir, distances)

  # find for each input arg if the internal point is inside or outside
  for (I,input_arg) in enumerate(input_args)
    bool_matrix[A,I]=is_internal_point(input_arg, internal_point, ray_dir)
  end

  if false
    points = GLBatch(POINTS); points.point_size = 8
    lines  = GLBatch(LINES) ; lines.line_width  = 2
    append!(points.colors.vector,RED); append!(points.vertices.vector, internal_point)
    for distance in distances
      append!(points.colors.vector, GREEN ); append!(points.vertices.vector, internal_point+distance*ray_dir)
      append!(lines.colors.vector,  YELLOW); append!(lines.vertices.vector, internal_point)
      append!(lines.colors.vector,  YELLOW); append!(lines.vertices.vector, internal_point+distance*ray_dir)
    end

    @show(join([string(it) for it in bool_matrix[A,:]], " "))
    text=GLText(join([string(it) for it in v], " "), center=internal_point, color=BLACK)
    VIEWCOMPLEX(atom, show=["V", "EV"], explode=[1.0,1.0,1.0], user_batches=[points,lines,text...])
  end
end

function Union(v)        return any(v) end
function Intersection(v) return all(v) end
function Difference(v)   return v[1] && !any(v[2:end]) end
function Xor(v)          return (length([it for it in v if it]) % 2)==1  end

bool_op=Intersection
sel=Cell()
for (A,row) in enumerate(eachrow(bool_matrix))
  if bool_op(row)
    append!(sel,atom_faces[A])
  end
end
sel=remove_duplicates(sel)
tmp=SELECT(arrangement, sel)
VIEWCOMPLEX(tmp, explode=[1.4,1.4,1.4])



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

# important to do simplification, because there are replicated faces
arrangement = SIMPLIFY(arrangement)

# show the atoms
show_atoms=false
if show_atoms
  @show(arrangement)
  for (A, sel) in enumerate(arrangement.C[:CF])
    atom = SELECT(arrangement, sel)
    @show(atom)
    VIEWCOMPLEX(atom, explode=[1.4, 1.4, 1.4], show=["V", "EV", "FV", "V_text", "EV_text", "FV_text"], face_color=TRANSPARENT)
  end
end

input_args=[LAR(it) for it in TOPOS(assembly)]
bool = bool3d(arrangement, bool_op=Difference, input_args=input_args, debug_mode=true)
VIEWCOMPLEX(bool, explode=[1.4,1.4,1.4])


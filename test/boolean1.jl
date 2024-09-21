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

arrangement=ARRANGE3D(lar)

atoms = ATOMS(arrangement, debug_mode=false)

input_args=[LAR(it) for it in TOPOS(assembly)]
bool = bool3d(atoms, bool_op=Difference, input_args=input_args, debug_mode=true)
VIEWCOMPLEX(bool, explode=[1.4,1.4,1.4])


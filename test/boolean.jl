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

input_args=[LAR(it) for it in TOPOS(assembly)]
bool = bool3d(arrangement, bool_op=Difference, input_args=input_args, debug_mode=true)
VIEWCOMPLEX(bool, explode=[1.4,1.4,1.4])


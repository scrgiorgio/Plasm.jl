using Plasm
using Random

Random.seed!(0)

cube = CUBOID([1,1,1])
hpc = STRUCT(
  cube, 
  T(1,2,3)(.5,.5,0.75), 
  R(2,3)(pi/4), 
  R(1,3)(pi/4),
  cube)
lar = LAR(hpc)

arrangement = ARRANGE3D(lar)
VIEWCOMPLEX(arrangement, explode=[1.2,1.2,2.0])
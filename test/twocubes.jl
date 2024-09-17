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

# show faces, exploding each face by its centroid
VIEWCOMPLEX(arrangement, show=["FV"], explode=[1.2,1.2,2.0])

# show faces, but keep the atom together
VIEWCOMPLEX(arrangement, show=["FV","atom"], explode=[1.2,1.2,2.0])
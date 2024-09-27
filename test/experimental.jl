using Plasm
using TetGen
using Random

Random.seed!(0)

hpc=STRUCT(CUBOID([1,1,1]), T(1,2,3)(0.5,0.5,0.5),CUBOID([1,1,1]))

hpc = STRUCT(
  CUBOID([1,1,1]), 
  RandomCube(0.2,2.0))

# hpc = STRUCT([RandomCube(0.2,2.0) for I in 1:6])

lar=LAR(hpc)

# arrangement=arrange3d(lar, debug_mode=false)

arrangement=fragment_lar(lar)

# show faces, exploding each face by its centroid
VIEWCOMPLEX(arrangement, show=["FV"], explode=[1.4,1.4,1.4])

# show faces, but keep the atom together
VIEWCOMPLEX(arrangement, show=["FV","atom"], explode=[1.4,1.4,1.4])

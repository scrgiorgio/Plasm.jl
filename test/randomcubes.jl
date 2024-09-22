using Plasm
using Random

Random.seed!(0)

hpc = STRUCT([RandomCube(0.2,2.0) for I in 1:6])
lar = LAR(hpc)

arrangement = ARRANGE3D(lar)

# remove outer
arrangement,___=SPLIT(arrangement) 

# show faces, exploding each face by its centroid
VIEWCOMPLEX(arrangement, show=["FV"], explode=[1.4,1.4,1.4])

# show faces, but keep the atom together
VIEWCOMPLEX(arrangement, show=["FV","atom"], explode=[1.4,1.4,1.4])
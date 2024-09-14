using Plasm
using Random

Random.seed!(0)

hpc = STRUCT([RandomCube(0.2,2.0) for I in 1:6])
lar = LAR(hpc)

arrangement = ARRANGE3D(lar)
VIEWCOMPLEX(arrangement; explode=[1.4,1.4,1.4])
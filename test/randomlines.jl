using Plasm
using Random

Random.seed!(0)

hpc = STRUCT([RandomLine(2.0,3.0) for I in 1:6])
lar = LAR(hpc)

arrangement = ARRANGE2D(lar)

# VIEWCOMPLEX(arrangement, show=["V","EV"], explode=[1.2,1.2,1.2])
VIEWCOMPLEX(arrangement, show=["V","FV"], explode=[1.2,1.2,1.2])
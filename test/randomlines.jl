using Plasm
using Random

Random.seed!(0)

hpc = STRUCT([RandomLine(2.0,3.0) for I in 1:6])
lar = LAR(hpc)

arrangement = ARRANGE2D(lar)

# show edges exploded by single edge
VIEWCOMPLEX(arrangement, show=["V","EV"], explode=[1.5,1.5,1.5])

# show edges exploded by atom centroid
VIEWCOMPLEX(arrangement, show=["V","EV", "atom"], explode=[1.5,1.5,1.5])

# show faces exploded by single face
VIEWCOMPLEX(arrangement, show=["V","FV"], explode=[1.5,1.5,1.5])

using Plasm
using Random

Random.seed!(0)

hpc = STRUCT([RandomBubble() for I in 1:50])
lar = LAR(hpc)

arrangement = ARRANGE2D(lar)

# show edges exploded by single edge
VIEWCOMPLEX(arrangement, show=["V","EV"], explode=[1.5,1.5,1.5])

# show edges, but explode to keep the face together
VIEWCOMPLEX(arrangement, show=["V","EV", "atom"], explode=[1.5,1.5,1.5])

# show faces exploded by single face
VIEWCOMPLEX(arrangement, show=["V","FV"], explode=[1.5,1.5,1.5])
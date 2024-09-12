using Plasm
using Random

Random.seed!(0)

hpc = STRUCT([RandomBubble() for I in 1:50])
lar = LAR(hpc)
V, EV  = lar.V, lar.C[:EV]
V,FVs,EVs = ARRANGE2D(V,EV)

VIEWCOMPLEX([Lar(V, Dict( :EV => EV)) for EV in EVs],explode=[1.2,1.2,1.2])
VIEWCOMPLEX([Lar(V, Dict( :FV => FV)) for FV in FVs],explode=[1.2,1.2,1.2])
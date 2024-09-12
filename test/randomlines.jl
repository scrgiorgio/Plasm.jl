using Plasm
using Random

Random.seed!(0)

hpc = STRUCT([RandomLine(2.0,3.0) for I in 1:6])
lar = LAR(hpc)
V,FVs,EVs = ARRANGE2D(lar.V,lar.C[:EV])

VIEWCOMPLEX([Lar(V, Dict( :EV => EV)) for EV in EVs], explode=[1.2,1.2,1.2])
VIEWCOMPLEX([Lar(V, Dict( :FV => FV)) for FV in FVs], explode=[1.2,1.2,1.2])
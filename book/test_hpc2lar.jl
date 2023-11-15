using Plasm
include("./hpc2lar.jl")


esa = CUBE(1)::Hpc
Lar(esa::Hpc)::Lar

twocubes = STRUCT([esa,T(3)(1),esa])::Hpc
doublecube = Lar(twocubes::Hpc)::Lar
KEV = lar2cop(doublecube.C[:EV])

update(doublecube::Lar, (:KEV, KEV))::Lar
doublecube.C
doublecube.V



using Plasm

# ////////////////////////////////////////////////////////////
function RandomCube(size_min::Float64,size_max::Float64)
  size = size_min+rand()*(size_max-size_min)
  return STRUCT(
    T(1,2,3)(rand(3)...), 
    S([1,2,3])([size,size,size]), 
    R([1,2])(2*pi*rand()),
    R([2,3])(2*pi*rand()),
    R([1,3])(2*pi*rand()),
    Plasm.CUBE(1) 
  )
end

# make sure I will get consistent random numbers (important for debugging)
import Random
Random.seed!(0)

hpc = STRUCT([RandomCube(0.2,2.0) for I in 1:6])
lar = LAR(hpc)

V, EV, FV  = lar.V, lar.C[:EV], lar.C[:FV]
V,CVs,FVs,EVs = testarrangement(V,FV,EV)
show_exploded(V,CVs,FVs,EVs)


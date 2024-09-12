using Plasm
using Random

Random.seed!(0)

cyl = T(3)(-2)(CYLINDER([0.5,4])(8))
hpc = STRUCT( 
  STRUCT(
    T(1,2,3)(-1,-1,-1),
    CUBE(2)), 
  cyl, R(2,3)(π/2), 
  cyl, R(1,3)(π/2), 
  cyl 
)
lar = LAR(hpc)

arrangement = ARRANGE3D(lar)
VIEWCOMPLEX(arrangement,explode=[1.2,1.2,1.8])
using Plasm
using Random

Random.seed!(0)

# /////////////////////////////////////////
function Run(primitive)

  # VIEW(primitive)

  hpc = STRUCT( 
    STRUCT(T(1,2,3)(-1,-1,-1),CUBE(2)), 
    primitive, R(2,3)(π/2), 
    primitive, R(1,3)(π/2), 
    primitive 
  )
  lar = LAR(hpc)
  # VIEWCOMPLEX(lar)
  
  arrangement = ARRANGE3D(lar)

  arrangement=without_outer_atom(arrangement)

  VIEWCOMPLEX(arrangement,explode=[1.2,1.2,1.8])
end

Run(T(3)(-2)(CYLINDER([0.5,4])(8)))
Run(T(3)(-2)(TUBE([0.3,0.5,4.0])(4)))


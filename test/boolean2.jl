using Plasm
using Random

Random.seed!(0)


primitive=T(3)(-2)(CYLINDER([0.5,4])(8))

hpc = STRUCT( 
  STRUCT(T(1,2,3)(-1,-1,-1),CUBE(2)), 
  primitive, R(2,3)(π/2), 
  primitive, R(1,3)(π/2), 
  primitive 
)
lar = LAR(hpc)

arrangement=ARRANGE3D(lar)

atoms=ATOMS(arrangement,debug_mode=false)

# any boolean expression will work
function my_bool_op(v)
  c,x1,x2,x3=v
  return Difference([c,Union([x1,x2,x3])])
end

boop_lar = bool3d(atoms, bool_op=my_bool_op, input_args=input_args, debug_mode=true)
VIEWCOMPLEX(boop_lar, explode=[1.4,1.4,1.4])

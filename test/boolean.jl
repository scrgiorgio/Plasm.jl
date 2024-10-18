using Plasm
using Random
using LinearAlgebra

Random.seed!(0)

# ////////////////////////////////////////////////////////
function TwoCubes()

  assembly = STRUCT(
    CUBE(1), 
    T(1,2,3)(.5,.5,.5), 
    CUBE(1)
  )

  lar=LAR(assembly)
  input_args=[LAR(it) for it in TOPOS(assembly)]
  lar=ARRANGE3D(lar)
  lar=INNERS(lar)
  lar= BOOL3D(lar, bool_op=Difference, input_args=input_args, debug_mode=false)
  @show(lar)
  VIEWCOMPLEX(lar, show=["FV"], explode=[1.4,1.4,1.4])
  # VIEWCOMPLEX(lar, show=["CV"], explode=[1.4,1.4,1.4]) # TODO inside boolean, SELECT does not produce CF

end

# ////////////////////////////////////////////////////////
function PieceCylinder()

  primitive=T(3)(-2)(CYLINDER([0.5,4])(8))
  assembly=  STRUCT( 
    STRUCT(T(1,2,3)(-1,-1,-1),CUBE(2)), 
    primitive, R(2,3)(π/2), 
    primitive, R(1,3)(π/2), 
    primitive 
  )

  # any boolean expression will work
  function my_bool_op(v)
    c,x1,x2,x3=v
    return Difference([c,Union([x1,x2,x3])])
  end

  lar = LAR(assembly)
  input_args=[LAR(it) for it in TOPOS(assembly)]
  lar=ARRANGE3D(lar)
  lar=INNERS(lar)

  VIEWCOMPLEX(lar)

  lar=BOOL3D(lar, bool_op=my_bool_op, input_args=input_args, debug_mode=false)
  VIEWCOMPLEX(lar, show=["FV"], explode=[1.4,1.4,1.4])
  # VIEWCOMPLEX(lar, show=["CV"], explode=[1.4,1.4,1.4]) # TODO inside boolean, SELECT does not produce CF

end

# ///////////////////////////////////////////
function Building()
  X = GRID([2.4,4.5,-3,4.5,2.4])
  Y = GRID([7,5])
  Z = GRID([3,3])
  assembly=  STRUCT(X * Y * Z)

  lar=LAR(assembly)
  input_args=[LAR(it) for it in TOPOS(assembly)]
  lar=ARRANGE3D(lar)
  lar=INNERS(lar)
  VIEWCOMPLEX(lar)
  lar= BOOL3D(lar, bool_op=Difference, input_args=input_args, debug_mode=false)
  @show(lar)
  VIEWCOMPLEX(lar, show=["FV"], explode=[1.4,1.4,1.4])
  # VIEWCOMPLEX(lar, show=["CV"], explode=[1.4,1.4,1.4]) # TODO inside boolean, SELECT does not produce CF

end

TwoCubes()
PieceCylinder()
# Building()

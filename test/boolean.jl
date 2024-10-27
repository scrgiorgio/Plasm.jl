using Plasm
using Random
using LinearAlgebra

Random.seed!(0)

# ////////////////////////////////////////////////////////
function TwoCubes()

  args=[
    CUBOID([0.0, 0.0, 0.0],[1.0, 1.0, 1.0]),
    CUBOID([0.5, 0.5, 0.5],[1.5, 1.5, 1.5])
  ]

  assembly=STRUCT(args)
  lar = ARRANGE3D(LAR(assembly))
  lar = INNERS(lar)

  for bool_op in [Union, Intersection, Difference, Xor]
    lar= BOOL3D(lar, bool_op=bool_op, input_args=[LAR(it) for it in args], debug_mode=false)
    VIEWCOMPLEX(lar, show=["FV"], explode=[1.0,1.0,1.0], title="TwoCubes FV $(bool_op)")
    VIEWCOMPLEX(lar, show=["CV"], explode=[1.2,1.2,1.2], title="TwoCubes CV $(bool_op)")
  end

end


# ////////////////////////////////////////////////////////
function ThreeCubes()

  args=[
    CUBOID([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
    CUBOID([0.5, 0.5, 0.0], [1.5, 1.5, 1.0]),
    CUBOID([0.7, 0.7, 0.0], [1.7, 1.7, 1.0])]

  lar=LAR(STRUCT(args))
  lar=INNERS(ARRANGE3D(lar))
  input_args=[LAR(arg) for arg in args]

  for bool_op in [Union, Intersection, Difference, Xor]
    lar=BOOL3D(lar, bool_op=bool_op, input_args=input_args)
    VIEWCOMPLEX(lar, show=["FV"], explode=[1.0,1.0,1.0], title="3cubes FV $(bool_op)")
    VIEWCOMPLEX(lar, show=["CV"], explode=[1.2,1.2,1.2], title="3cubes CV $(bool_op)")
  end
  
end

# ////////////////////////////////////////////////////////
function PieceCylinder()

  cyl=T(3)(-2)(CYLINDER([0.5,4])(8))
  cube=CUBOID([-1,-1,-1],[+1,+1,+1])

  args=[
    cube,
    R(1,2)(π/2)(cyl),
    R(2,3)(π/2)(cyl),
    R(1,3)(π/2)(cyl)
  ]

  assembly=  STRUCT(args)

  # any boolean expression will work
  #function my_bool_op(v)
  #  c,x1,x2,x3=v
  #  return Union([c,Union([x1,x2,x3])])
  #end

  lar = LAR(assembly)
  input_args=[LAR(it) for it in TOPOS(assembly)]
  lar=ARRANGE3D(lar)
  lar=INNERS(lar)

  for bool_op in [Union, Intersection, Difference, Xor]
    lar=BOOL3D(lar, bool_op=my_bool_op, input_args=input_args, debug_mode=false)
    VIEWCOMPLEX(lar, show=["FV"], explode=[1.0,1.0,1.0], title="PieceCylinder FV")
    VIEWCOMPLEX(lar, show=["CV"], explode=[1.2,1.2,1.2], title="PieceCylinder CV")
  end

end

TwoCubes()
# ThreeCubes()
# PieceCylinder()

# does not make sense because some cells return "outside" , some other "inside depending or where I am inside the grid
# Building()

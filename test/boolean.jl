using Plasm
using Random
using LinearAlgebra

Random.seed!(0)

# /////////////////////////////////////////////////////////////////////////////
function RunBooleanTest(name, bool_op, args; debug_mode=false)
  assembly=STRUCT(args)
  lar = ARRANGE3D(LAR(assembly))
  lar = INNERS(lar)

  if debug_mode
    for atom in ATOMS(lar)
      show_debug(atom)
      VIEWCOMPLEX(atom, explode=[1.0, 1.0, 1.0], show=["V", "EV","FV"])
    end
  end

  result= BOOL3D(lar, bool_op=bool_op, input_args=[LAR(arg) for arg in args], debug_mode=debug_mode)
  VIEWCOMPLEX(result, show=["FV"], explode=[1.2,1.2,1.2], title="$(name)/$(bool_op) FV")
  VIEWCOMPLEX(result, show=["CV"], explode=[1.2,1.2,1.2], title="$(name)/$(bool_op) CV")
end

# //////////////////////////////////////////////////////
function TwoCubes()
  return [
    CUBOID([0.0, 0.0, 0.0],[1.0, 1.0, 1.0]),
    CUBOID([0.5, 0.5, 0.5],[1.5, 1.5, 1.5])
  ]
end

# //////////////////////////////////////////////////////
function ThreeCubes()
  return [
    CUBOID([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
    CUBOID([0.5, 0.5, 0.0], [1.5, 1.5, 1.0]),
    CUBOID([0.7, 0.7, 0.0], [1.7, 1.7, 1.0])
    ]
end

# //////////////////////////////////////////////////////
function MechanicalPiece(primitive::Hpc)
  return [
    CUBOID(
      [-1.0,-1.0,-1.0],
      [+1.0,+1.0,+1.0]
    ),
    R(1,2)(π/2)(primitive),
    R(2,3)(π/2)(primitive),
    R(1,3)(π/2)(primitive)
  ]
end

# //////////////////////////////////////////////////////
function MechanicalPiece(symbol::Symbol)
  if symbol==:cube
    return MechanicalPiece(CUBOID( [-0.4, -0.4, -2.0], [+0.4, +0.4, +2.0]))
  elseif symbol==:cylinder
    return MechanicalPiece(T(3)(-2)(CYLINDER([0.5,4])(8)))
  elseif symbol==:tube
    return MechanicalPiece(T(3)(-2)(TUBE([0.3,0.5,4.0])(4)))
  else
    @assert(false)
  end
end


# //////////////////////////////////////////////////////
#for bool_op in [Union, Intersection, Difference, Xor]
for bool_op in [Union, Intersection, Difference, Xor]
  RunBooleanTest("boolean/TwoCubes",bool_op, TwoCubes())
end

for bool_op in [Union, Intersection, Difference, Xor]
  RunBooleanTest("boolean/ThreeCubes",bool_op, ThreeCubes())
end

# example of user-defined boolean expression: any boolean expression will work
function UserBoolOp(v::Vector{Bool})::Bool
  return Difference([v[1],Union(v[2:end])])
end

for bool_op in [Union, Intersection, Difference, Xor]
  RunBooleanTest("boolean/MechanicalPiece/cube", bool_op, MechanicalPiece(:cube))
end

for bool_op in [Union, Intersection, Difference, Xor]
  RunBooleanTest("boolean/MechanicalPiece/cylinder", bool_op, MechanicalPiece(:cylinder))
end

# BROKEN
if false
  for bool_op in [Union, Intersection, Difference, Xor]
    RunBooleanTest("boolean/MechanicalPiece/tube", bool_op, MechanicalPiece(:tube))
  end
end

# does not make sense because some cells return "outside" , some other "inside depending or where I am inside the grid
# Building()

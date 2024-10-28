using Plasm
using Random


# ////////////////////////////////////////////
function TwoCubes()
  return STRUCT(
    CUBOID([1,1,1]), 
    T(1,2,3)(0.5,0.5,0.5),
    CUBOID([1,1,1]))
end

# ////////////////////////////////////////////
function RandomCubes(ncubes=6)
  return STRUCT([RandomCube(0.2,2.0) for I in 1:ncubes])
end

# ////////////////////////////////////////////
function MechanicalPieceCylinder(num_subdivisions::Int=8)
  primitive=T(3)(-2)(CYLINDER([0.5,4])(num_subdivisions))
  return STRUCT( 
    STRUCT(T(1,2,3)(-1,-1,-1),CUBE(2)), 
    primitive, R(2,3)(π/2), 
    primitive, R(1,3)(π/2), 
    primitive 
  )
end

# ////////////////////////////////////////////
function MechanicalPieceTube(num_subdivisions::Int=4)
  primitive=T(3)(-2)(TUBE([0.3,0.5,4.0])(num_subdivisions))
  return STRUCT( 
    STRUCT(T(1,2,3)(-1,-1,-1),CUBE(2)), 
    primitive, R(2,3)(π/2), 
    primitive, R(1,3)(π/2), 
    primitive 
  )
end

# ///////////////////////////////////////////
function Building()
  X = GRID([2.4,4.5,-3,4.5,2.4])
  Y = GRID([7,5])
  Z = GRID([3,3])
  return  X * Y * Z
end

# ///////////////////////////////////////////
function RunArrange3DTest(name, hpc::Hpc)
  lar=LAR(hpc)
  lar=ARRANGE3D(lar)
  lar=INNERS(lar)
  VIEWCOMPLEX(lar, show=["FV"], explode=[1.4,1.4,1.4], title="$(name) FV")
  VIEWCOMPLEX(lar, show=["CV"], explode=[1.4,1.4,1.4], title="$(name) CV")
end

# //////////////////////////////////////////////
begin
  Random.seed!(0)

  RunArrange3DTest("arrange3d/twocubes", TwoCubes())
  RunArrange3DTest("arrange3d/building", Building())
  RunArrange3DTest("arrange3d/randomcubes", RandomCubes())
  RunArrange3DTest("arrange3d/MechanicalPiece/cylinder", MechanicalPieceCylinder())
  RunArrange3DTest("arrange3d/MechanicalPiece/tube", MechanicalPieceTube())

end


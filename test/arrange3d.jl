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
function PieceCylinder(num_subdivisions::Int=8)
  primitive=T(3)(-2)(CYLINDER([0.5,4])(num_subdivisions))
  return STRUCT( 
    STRUCT(T(1,2,3)(-1,-1,-1),CUBE(2)), 
    primitive, R(2,3)(π/2), 
    primitive, R(1,3)(π/2), 
    primitive 
  )
end

# ////////////////////////////////////////////
function PieceTube(num_subdivisions::Int=4)
  primitive=T(3)(-2)(TUBE([0.3,0.5,4.0])(num_subdivisions))
  return STRUCT( 
    STRUCT(T(1,2,3)(-1,-1,-1),CUBE(2)), 
    primitive, R(2,3)(π/2), 
    primitive, R(1,3)(π/2), 
    primitive 
  )
end

# ///////////////////////////////////////////
function View3D(hpc::Hpc)
  lar=LAR(hpc)
  lar=ARRANGE3D(lar)
  lar=INNERS(lar)
  VIEWCOMPLEX(lar, show=["FV"], explode=[1.4,1.4,1.4], title="Inner atoms")
  VIEWCOMPLEX(lar, show=["FV","atom"], explode=[1.4,1.4,1.4], title="Inner atoms")
end

# //////////////////////////////////////////////
begin
  Random.seed!(0)

  #Plasm.LAR_ARRANGE_VERSION=1
  #View3D(TwoCubes())
  #View3D(RandomCubes())
  #View3D(PieceCylinder())
  #View3D(PieceTube())

  Plasm.LAR_ARRANGE_VERSION=2
  hpc=TwoCubes()
  # hpc=RandomCubes(2)
  # hpc=RandomCubes(6)
  # hpc=PieceCylinder()
  # hpc=PieceTube()
  lar=LAR(hpc)
  lar=ARRANGE3D(lar, debug_mode=true)
  lar=INNERS(lar)
  VIEWCOMPLEX(lar, show=["FV","atom"], explode=[1.2,1.2,1.2])

end


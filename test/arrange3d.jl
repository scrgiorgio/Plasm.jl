using Plasm
using TetGen
using Random

Random.seed!(0)

function TwoCubes()
  return STRUCT(
    CUBOID([1,1,1]), 
    T(1,2,3)(0.5,0.5,0.5),
    CUBOID([1,1,1]))
end


function RandomCubes(ncubes=6)
  return STRUCT([RandomCube(0.2,2.0) for I in 1:ncubes])
end

function PieceCylinder(num_subdivisions::Int=8)
  primitive=T(3)(-2)(CYLINDER([0.5,4])(num_subdivisions))
  return STRUCT( 
    STRUCT(T(1,2,3)(-1,-1,-1),CUBE(2)), 
    primitive, R(2,3)(π/2), 
    primitive, R(1,3)(π/2), 
    primitive 
  )
end

function PieceTube(num_subdivisions::Int=4)
  primitive=T(3)(-2)(TUBE([0.3,0.5,4.0])(num_subdivisions))
  return STRUCT( 
    STRUCT(T(1,2,3)(-1,-1,-1),CUBE(2)), 
    primitive, R(2,3)(π/2), 
    primitive, R(1,3)(π/2), 
    primitive 
  )
end

function View3D(lar::Lar)
  VIEWCOMPLEX(lar, show=["FV"], explode=[1.4,1.4,1.4])
  VIEWCOMPLEX(lar, show=["FV","atom"], explode=[1.4,1.4,1.4])
end

View3D(SPLIT(ARRANGE3D(LAR(TwoCubes()     )))[1])
View3D(SPLIT(ARRANGE3D(LAR(RandomCubes()  )))[1])
View3D(SPLIT(ARRANGE3D(LAR(PieceCylinder())))[1])
View3D(SPLIT(ARRANGE3D(LAR(PieceTube()    )))[1])




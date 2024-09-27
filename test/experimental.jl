using Plasm
using Random

Random.seed!(0)

# ///////////////////////////////////
function RandomSquares(num_squares::Int=6)
  return STRUCT([RandomSquare(2.0,3.0) for I in 1:num_squares ])
end

# ///////////////////////////////////
function RandomBubbles(num_bubbles::Int=50)
  return STRUCT([RandomBubble() for I in 1:num_bubbles])
end

# ///////////////////////////////////
function TwoCubes()
  return STRUCT(
    CUBOID([1,1,1]), 
    T(1,2,3)(0.5,0.5,0.5),
    CUBOID([1,1,1]))
end

# ///////////////////////////////////
function RandomCubes(ncubes=6)
  return STRUCT([RandomCube(0.2,2.0) for I in 1:ncubes])
end

# ///////////////////////////////////
function PieceCylinder(num_subdivisions::Int=8)
  primitive=T(3)(-2)(CYLINDER([0.5,4])(num_subdivisions))
  return STRUCT( 
    STRUCT(T(1,2,3)(-1,-1,-1),CUBE(2)), 
    primitive, R(2,3)(π/2), 
    primitive, R(1,3)(π/2), 
    primitive 
  )
end

# ///////////////////////////////////
function PieceTube(num_subdivisions::Int=4)
  primitive=T(3)(-2)(TUBE([0.3,0.5,4.0])(num_subdivisions))
  return STRUCT( 
    STRUCT(T(1,2,3)(-1,-1,-1),CUBE(2)), 
    primitive, R(2,3)(π/2), 
    primitive, R(1,3)(π/2), 
    primitive 
  )
end

# ///////////////////////////////////
function View2D(lar::Lar)
  VIEWCOMPLEX(lar, show=["V","EV"        ], explode=[1.5,1.5,1.5])
  VIEWCOMPLEX(lar, show=["V","EV", "atom"], explode=[1.5,1.5,1.5])
  VIEWCOMPLEX(lar, show=["V","FV"        ], explode=[1.5,1.5,1.5])
end

# ///////////////////////////////////
function View3D(lar::Lar)
  VIEWCOMPLEX(lar, show=["FV"       ], explode=[1.4,1.4,1.4])
  VIEWCOMPLEX(lar, show=["FV","atom"], explode=[1.4,1.4,1.4])
end

# View2D(arrange2d_experimental(LAR(RandomSquares())))
# View2D(arrange2d_experimental(LAR(RandomBubbles())))

# View3D(fragment_lar(LAR(TwoCubes())))
# View3D(fragment_lar(LAR(RandomCubes())))
# View3D(fragment_lar(LAR(PieceCylinder())))
# View3D(fragment_lar(LAR(PieceTube())))

# broken
# View3D(arrange3d_experimental(LAR(TwoCubes())))
# View3D(arrange3d_experimental(LAR(RandomCubes())))
# View3D(arrange3d_experimental(LAR(PieceCylinder())))
# View3D(arrange3d_experimental(LAR(PieceTube())))

# arrange3d(lar, debug_mode=true)



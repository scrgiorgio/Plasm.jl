using Plasm
using Random

# ////////////////////////////////////////////////
function RandomSquares(num_squares::Int=6)
  return STRUCT([RandomSquare(2.0,3.0) for I in 1:num_squares ])
end

# ////////////////////////////////////////////////
function RandomBubbles(num_bubbles::Int=50)
  return STRUCT([RandomBubble() for I in 1:num_bubbles])
end

# ////////////////////////////////////////////////
function View2D(lar::Lar; explode=[1.0,1.0,1.0])
  #VIEWCOMPLEX(lar, show=["V","EV"],         explode=[1.5,1.5,1.5])
  #VIEWCOMPLEX(lar, show=["V","EV", "atom"], explode=[1.5,1.5,1.5])
  VIEWCOMPLEX(lar, show=["V","FV"],          explode=explode)
end

# ////////////////////////////////////////////////
begin

  Random.seed!(0)

  Plasm.LAR_ARRANGE_VERSION=1
  #View2D(ARRANGE2D(LAR(RandomSquares())))
  #View2D(ARRANGE2D(LAR(RandomBubbles())))

  Plasm.LAR_ARRANGE_VERSION=2
  View2D(ARRANGE2D(LAR(RandomSquares())))
  View2D(ARRANGE2D(LAR(RandomBubbles())))
  
end

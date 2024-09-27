using Plasm
using Random

Random.seed!(0)

function RandomSquares(num_squares::Int=6)
  return STRUCT([RandomSquare(2.0,3.0) for I in 1:num_squares ])
end

function RandomBubbles(num_bubbles::Int=50)
  return STRUCT([RandomBubble() for I in 1:num_bubbles])
end

function View2D(lar::Lar)
  VIEWCOMPLEX(lar, show=["V","EV"], explode=[1.5,1.5,1.5])
  VIEWCOMPLEX(lar, show=["V","EV", "atom"], explode=[1.5,1.5,1.5])
  VIEWCOMPLEX(lar, show=["V","FV"], explode=[1.5,1.5,1.5])
end

View2D(ARRANGE2D(LAR(RandomSquares())))
View2D(ARRANGE2D(LAR(RandomBubbles())))


using Plasm
using Random
using LinearAlgebra

Random.seed!(0)

SK = SKELETON

X = GRID([2.4,4.5,-3,4.5,2.4])
Y = GRID([7,5])
Z = GRID([3,3])

idea = X * Y * Z

assessement = (
  T(3)(.5) ∘ 
  S(1,2,3)(1+1/91,1+1/60,-1)
)

grounds =  assessement(
  BOX([1, 2])(idea)* 
  GRID([0.2, -6, 1])
)

build1_110=(SK(1)(X*Y)*SK(0)(Z))
build1_101=(SK(0)(X)*SK(1)(Y*Z))
build1_011=(R(2,3)(-π/2)(SK(1)(X*Z)*SK(0)(Y)))

structural_frames = [
  build1_110, 
  build1_101, 
  build1_011]

hpc = STRUCT(
  AA(OFFSET([.2,.2,.5]))(structural_frames), 
  grounds)

VIEW(hpc)

# BROKEN (!)
#arrangement=ARRANGE3D(LAR(hpc))
#arrangement=INNERS(arrangement)
#VIEW(outer)
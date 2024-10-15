using Plasm
SK = SKELETON

X = GRID([2.4,4.5,-3,4.5,2.4])
Y = GRID([7,5])
Z = GRID([3,3])
idea = X * Y * Z

building110 = X*Y*SK(0)(Z)
building101 = X*SK(0)(Y)*Z
building011 = SK(0)(X)*Y*Z

building1_110 = SK(1)(X*Y)*SK(0)(Z)
building1_101 = SK(0)(X)*SK(1)(Y*Z)
building1_011 = SK(1)(X*Z)*SK(0)(Y)

dy = SIZE(2)(building1_011)    
building1_011 = (R(2,3)(-π/2))(building1_011)
building1_011 = (T(3)(dy))(building1_011)

floors = OFFSET([.2,.2,.2])(building110)
framex = OFFSET([.2,.2,.2])(building1_011)
framey = OFFSET([.2,.2,-.4])(building1_101)
framexy = STRUCT(framex, framey)
framexyz = STRUCT(framex, framey, floors)
VIEWCOMPLEX(LAR(framexyz))

VIEWCOMPLEX(LAR(SK(1)(idea)))




using Plasm

X = GRID([2.4,4.5,-3,4.5,2.4]);
Y = GRID([7,5]); Z = GRID([3,3]); 
idea = X * Y * Z;
#VIEWCOMPLEX( LAR(idea), explode = [1.2,1.2,2.0] )

building110 = X*Y*SK(0)(Z); 
#VIEWCOMPLEX(LAR(building110))  
building101 = X*SK(0)(Y)*Z; 
#VIEWCOMPLEX(LAR(building101)) 
building011 = SK(0)(X)*Y*Z; 
#VIEWCOMPLEX(LAR(building011)) 

building1_110 = SK(1)(X*Y)*SK(0)(Z);
#VIEWCOMPLEX(LAR(building1_110));  VIEW(building1_110)
building1_101 = SK(0)(X)*SK(1)(Y*Z);
#VIEWCOMPLEX(LAR(building1_101)); VIEW(building1_101)
building1_011 = SK(1)(X*Z)*SK(0)(Y); 
#VIEW(building1_011)
dy = SIZE(2)(building1_011);        
building1_011 = (R(2,3)(-π/2))(building1_011)
#VIEW(building1_011)
building1_011 = (T(3)(dy))(building1_011)
#VIEW(building1_011)

floors = OFFSET([.2,.2,.2])(building110);
#VIEWCOMPLEX(LAR(floors))
framex = OFFSET([.2,.2,.2])(building1_011);
#VIEWCOMPLEX(LAR(framex))
framey = OFFSET([.2,.2,-.4])(building1_101);
#VIEWCOMPLEX(LAR(framey))
framexy = STRUCT(framex, framey);
#VIEWCOMPLEX(LAR(framexy))
framexyz = STRUCT(framex, framey, floors);
VIEWCOMPLEX(LAR(framexyz))
#arrangement = ARRANGE3D(LAR(framexyz))
#VIEWCOMPLEX(arrangement, show("FV","atom"))



using Plasm

# //////////////////////////////////////////////////////////////////////////////
# 3D Boolean example generation (see CAD23 paper)
cube = CUBE(1)
assembly = STRUCT(cube, 
   R(1,3)(pi/4), T(1,2,3)(.5,.5,.5), cube, 
   R(2,3)(pi/5), T(1,2,3)(.5,.5,.5), cube)
W, (copEV, copFE, copCF), boolmatrix = bool3d(assembly);

V,CVs,FVs,EVs = pols2tria(W, copEV, copFE, copCF);
show_exploded(V,CVs,FVs,EVs)

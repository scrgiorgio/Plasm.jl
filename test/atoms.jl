using Plasm, SparseArrays, DataStructures, LinearAlgebra
# generate a 3D space arrangement
cube = CUBE(1)
assembly = STRUCT(cube, R(1,2)(π/4), cube)	
global V,FV,EV = new2old(LAR(assembly)) # V by columns
# initialization
cop_EV = convert(ChainOp, coboundary_0(EV::Cells))
cop_FE = coboundary_1(V, FV::Cells, EV::Cells); ## TODO: debug
W = convert(Points, V')
# generate the 3D space arrangement
V, copEV, copFE, copCF = space_arrangement( W, cop_EV, cop_FE )
@show (V, copEV, copFE, copCF);
# generate and draw the atoms [and the b-rep of outer space]
V,pols,_ = chainbasis2solids(V,copEV,copFE,copCF)
#@show V,pols;
DRAWATOMS(V,copEV,copFE,copCF, pols;outer=true)

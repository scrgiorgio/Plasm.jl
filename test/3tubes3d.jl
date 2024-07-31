using Plasm

tube = T(3)(-2)(TUBE([0.5,0.4,4])(4));
assembly = STRUCT( tube, R(2,3)(π/2), tube, R(1,3)(π/2), tube );
VIEW(assembly)

obj = LAR(assembly);
V, (EV, FV) = obj.V, (obj.C[:EV], obj.C[:FV]);

W = permutedims(V);
KEV = lar2cop(EV);
KFV = lar2cop(FV);

KFE = (KFV * KEV') .÷ Int8(2);
V, copEV, copFE, copCF = space_arrangement(
   W::Points, KEV::ChainOp, KFE::ChainOp);
   
# va in loop durante lo splitting 2D, scrivendo:
# 1/60
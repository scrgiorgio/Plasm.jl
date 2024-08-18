using Plasm  # da usare per testare TGW:  sembra terminare correttamente,      ma senza la condizione sufficiente "all(marks[k]=2 for all k"

rod = T(3)(-2)(CYLINDER([0.5,4])(4));
assembly = STRUCT( rod, R(2,3)(π/2), rod, R(1,3)(π/2), rod );
VIEW(assembly)

obj = LAR(assembly);
V, (EV, FV) = obj.V, (obj.C[:EV], obj.C[:FV]);

V,CVs,FVs,EVs = testarrangement(V,FV,EV)
show_exploded(V,CVs,FVs,EVs)


#tube = T(3)(-2)(TUBE([0.5,0.4,4])(4));

#W = permutedims(V);
#KEV = lar2cop(EV);
#KFV = lar2cop(FV);
#
#KFE = (KFV * KEV') .÷ Int8(2);
#V, copEV, copFE, copCF = space_arrangement(
#   W::Points, KEV::ChainOp, KFE::ChainOp);
#   
## va in loop durante lo splitting 2D, scrivendo:
## 1/60
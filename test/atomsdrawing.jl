
using Plasm, SparseArrays
cube = CUBE(1)
assembly = STRUCT(cube, R(1,2)(π/4), cube)	

   ensemble = LAR(assembly)
   V,FV,EV = new2old(ensemble)
   #----------------------------------------------------------------------------
	#V,FV,EV = struct2lar(assembly)
	cop_EV = convert(ChainOp, coboundary_0(EV::Cells));
	cop_FE = coboundary_1(V, FV::Cells, EV::Cells); ## TODO: debug
	W = convert(Points, V');
	# generate the 3D space arrangement
	#----------------------------------------------------------------------------
	# generate the 3D space arrangement
	V, copEV, copFE, copCF = space_arrangement( W, cop_EV, cop_FE );
	@show V, copEV, copFE, copCF;


EV = AA(sort)([findnz(copEV[k,:])[1] for k=1:copEV.m]) # vertices per edge
FE = AA(sort)([findnz(copFE[k,:])[1] for k=1:copFE.m]) # edges per face
FV = AA(sort)([union(CAT([EV[e] for e in f])) for f in FE]) # vertices per face

dictEV = Dict{Array{Int},Int}(collect(zip(EV,1:length(EV))))
dictFV = Dict{Array{Int},Int}(collect(zip(FV,1:length(FV))))

# Visualization of whole arrangement
W = [V[:,k] for k=1:size(V,2)]
VIEWCOMPLEX(LAR(MKPOL(W,FV,EV))) # Intero modello dopo l'arrangement


# extraction of atoms from arrangement
U,pols,CF = chainbasis2solids(V,copEV,copFE,copCF)

# Visualization of all atoms
for k=1:length(pols)
   EV,FV = pols[k]; EV=CAT(EV)
   VIEWCOMPLEX(LAR(MKPOL(U,FV,EV)))
end

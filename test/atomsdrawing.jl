
using Plasm, SparseArrays, DataStructures

cube = CUBE(1)
assembly = STRUCT(cube, R(1,2)(π/4), cube)	

ensemble = LAR(assembly)
global V,FV,EV = new2old(ensemble)

#----------------------------------------------------------------------------
#V,FV,EV = struct2lar(assembly)
cop_EV = convert(ChainOp, coboundary_0(EV::Cells));
cop_FE = coboundary_1(V, FV::Cells, EV::Cells); ## TODO: debug
W = convert(Points, V');
#----------------------------------------------------------------------------

# generate the 3D space arrangement
V, copEV, copFE, copCF = space_arrangement( W, cop_EV, cop_FE );
#@show V, copEV, copFE, copCF;

# Lar complex computation, with associated dictionaries
EV = AA(sort)([findnz(copEV[k,:])[1] for k=1:copEV.m]) # vertices per edge
FE = AA(sort)([findnz(copFE[k,:])[1] for k=1:copFE.m]) # edges per face
FV = AA(sort)([union(CAT([EV[e] for e in f])) for f in FE]) # vertices per face
W = [V[k,:] for k=1:size(V,1)]

Wdict = OrderedDict{Any,Int}(zip(W,1:LEN(W)));
dictEV = Dict{Vector{Int},Int}(collect(zip(EV,1:length(EV))))
dictFV = Dict{Vector{Int},Int}(collect(zip(FV,1:length(FV))))

# Visualization of whole arrangement
arranged = LAR(MKPOL(W,FV,EV))

# Intero modello dopo l'arrangement
VIEWCOMPLEX(arranged)  

# extraction of atoms from arranged space
_,pols,_ = chainbasis2solids(V,copEV,copFE,copCF)

# Visualization of all atoms
for k=1:length(pols)

   localV=permutedims(V)
   localEV,localFV = pols[k] 
   localEV = union(CAT(localEV))
   localFV = AA(sort)(localFV)
   
   VIEWCOMPLEX2(
      localV, localEV, localFV, 
      [string(Wdict[a])  for a in W],
      [string(dictEV[a]) for a in localEV],
      [string(dictFV[a]) for a in localFV],
      properties=Properties("background_color"=>DEFAULT_BACKGROUND_COLOR))
  
end

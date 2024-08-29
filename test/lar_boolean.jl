using Plasm, SparseArrays

# //////////////////////////////////////////////////////////////////////////////
# 3D Boolean example generation (see CAD23 paper)
cube = CUBE(1);
assembly = STRUCT(cube, 
   R(1,3)(pi/4), T(1,2,3)(.5,.5,.5), cube, 
   R(2,3)(pi/5), T(1,2,3)(.5,.5,.5), cube);
#VIEW(assembly);
   V,FV,EV = new2old(LAR(assembly)) # shape(V) = 3×24 Matrix{Float64}:
   #----------------------------------------------------------------------------
	#V,FV,EV = struct2lar(assembly)
	cop_EV = convert(ChainOp, coboundary_0(EV::Cells));
	cop_FE = coboundary_1(V, FV::Cells, EV::Cells); ## TODO: debug
	W = convert(Points, V');
	# generate the 3D space arrangement
	#----------------------------------------------------------------------------
	# generate the 3D space arrangement
	V, copEV, copFE, copCF = space_arrangement(W, cop_EV, cop_FE );
	
boolmatrix = bool3d(assembly, V,copEV,copFE,copCF);







V,CVs,FVs,EVs = pols2tria(W, copEV, copFE, copCF);
show_exploded(V,CVs,FVs,EVs)

Matrix(boolmatrix)
#three-chains = [ for k = 1:3]

A = boolmatrix[:,2];
B = boolmatrix[:,3];
C = boolmatrix[:,4];
AorB = A .| B;
AandB = A .& B;
AxorB = AorB .⊻ (.!AandB) # = A .⊻ B;
AorBorC = A .| B .| C
AorBorC = .|(A, B, C)
AandBandC = A .& B .& C
AandBandC = .&(A, B, C)
AminBminC = .&(A, .!B, .!C) # A - B - C

unione = Matrix(copCF)' * Int.(AorBorC); # coord vector of Faces
intersection = Matrix(copCF)' * Int.(AandBandC); # coord vector of Faces
difference = Matrix(copCF)' * Int.(AminBminC); # coord vector of Faces

V,CVs,FVs,EVs = pols2tria(W, copEV, copFE, copCF) # whole assembly
Fs = difference
V,CVs,FVs,EVs = pols2tria(W, copEV, copFE, copCF, Fs) # part of assembly
V,CVs,FVs,EVs = pols2tria(W, copEV, copFE, copCF, intersection) # subassembly
V,CVs,FVs,EVs = pols2tria(W, copEV, copFE, copCF, unione) # subassembly


#EV = cop2lar(copEV)
#FE = cop2lar(copFE)
#CF = cop2lar(convert(ChainOp,copCF))
#EVor = [ev for (k,ev) in enumerate(EV) if abs(unione[k])==1 ]
#EVand = [ev for (k,ev) in enumerate(EV) if abs(intersection[k])==1 ]
#EVxor = [ev for (k,ev) in enumerate(EV) if abs(xor[k])==1 ]

using ViewerGL
GL = ViewerGL
GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1));
GL.VIEW(GL.GLExplode(V,EVs,1.,1.,1.,1,1));
#meshes = GL.GLExplode(V,CVs[1:end],5,5,5,99,1);
#GL.VIEW( push!( meshes, GL.GLFrame) );

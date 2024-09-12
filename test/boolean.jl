using Plasm
using Random

Random.seed!(0)

assembly = STRUCT(
  CUBE(1), 
  T(1,2,3)(.5,.5,.5), 
  CUBE(1)
)

lar=LAR(assembly)

arrangement = ARRANGE3D(lar)
VIEWCOMPLEX(arrangement,explode=[1.4,1.4,1.4])

bool3d(assembly, arrangement, CF)

"""

# 3 atoms
A = boolmatrix[:,2];
B = boolmatrix[:,3];
C = boolmatrix[:,4];
AorB = A .| B;
AandB = A .& B;
AxorB = AorB .⊻ (.!AandB) 
AorBorC = A .| B .| C
AorBorC = .|(A, B, C)
AandBandC = A .& B .& C
AandBandC = .&(A, B, C)
AminBminC = .&(A, .!B, .!C) # A - B - C

union  = Matrix(copCF)' * Int.(AorBorC  ) # coord vector of Faces
inters = Matrix(copCF)' * Int.(AandBandC) # coord vector of Faces
diff   = Matrix(copCF)' * Int.(AminBminC) # coord vector of Faces

arrangement = VIEWCOMPLEX(V_original, copEV, copFE, copCF) 
arrangement = VIEWCOMPLEX((SELECT(V_original, copEV, copFE, copCF, diff  )...)
arrangement = VIEWCOMPLEX(SELECT(V_original, copEV, copFE, copCF, inters) ...)
arrangement = VIEWCOMPLEX((SELECT(V_original, copEV, copFE, copCF, union )...)

"""
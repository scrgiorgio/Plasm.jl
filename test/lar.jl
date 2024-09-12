using Plasm

import Random

# ///////////////////////////////////////////////////////
function TestToLAR()

  Random.seed!(0)

  lar = LAR(SPHERE(1.)([3,3]));
  VIEW(MKPOLS(lar.V,lar.C[:EV]))
  VIEW(MKPOLS(lar.V,lar.C[:FV]))


	# test2D
  begin
    points=[
      [0.5,0.5],
      [0.0,0.0],[1.0,0.0],[1.0,1.0], [0.5,2.0], [0.0,1.0],
      [0.6,0.6],
      [0.0,0.0],[1.0,0.0],[1.0,1.0], [0.5,2.0], [0.0,1.0],
      [0.9,0.2]
    ]
    hpc=MkPol(points)
    geo = ToGeometry(hpc)
    @assert length(geo.points)==5 
    @assert length(geo.hulls)==1
    @assert length(geo.faces)==1
  end

	# test3D
  begin
    points=[
      [0.5,0.5,0.5],
      [0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0], [0.5,2.0,0.0], [0.0,1.0,0.0],
      [0.6,0.6,0.6],
      [0.0,0.0,1.0],[1.0,0.0,1.0],[1.0,1.0,1.0], [0.5,2.0,1.0], [0.0,1.0,1.0],
      [0.9,0.2,0.3]
    ]
    hpc=MkPol(points)
    geo = ToGeometry(hpc)
    @assert length(geo.points)==10 
    @assert length(geo.hulls)==1
    @assert length(geo.faces)==7
  end

  # test two 3d cubes
  begin
    hpc=STRUCT([
      CUBOID([1,1,1]),
      T([1])([3]),
      CUBOID([1,1,1])
    ])

    geo = ToGeometry(hpc)
    @assert(length(geo.points)==8*2)    # each cube is 8 vertices
    @assert(length(geo.hulls )==2)       # each cube is 1 hull
    @assert(length(geo.faces)==6*2)    # each cube has 6 boundary faces
  end

  # 2D b-rep
  begin                                   
    geo=ToGeometry(CIRCUMFERENCE(1.0)(8)) 
    @assert(length(geo.points[1])==2)     
    @assert(length(geo.hulls)==0)         
    @assert(length(geo.edges)==8)         
    @assert(length(geo.hulls)==0)         
  end                                     
  

  # 3D brep
  begin                                          
    geo=ToGeometry(SPHERE(1.0)([4,8]),precision=0.0)
    @assert(length(geo.points)==(4+1)*(8+1))     
    @assert(length(geo.points[1])==3)            
    @assert(length(geo.faces)==4*8)              
    @assert(length(geo.hulls)==0)                
  end   
  
  
  # show how to keep all cells in Hpc/Lar transformation
  begin  
    hpc1=STRUCT(
      T(1)(0.1),T(2)(0.1),T(3)(0.1),
      T(1)(0.0), Simplex(1),
      T(1)(1.1), Simplex(2),
      T(1)(1.1), Simplex(3),
      T(1)(1.1), Cube(1),
      T(1)(1.1), Cube(2),
      T(1)(1.1), Cube(3),
      T(1)(2.1), CIRCUMFERENCE(1.0)(8),   
      T(1)(2.1), SPHERE(1.0)([4,8])
    )            

    lar=LAR(hpc1)

    # to keep edges I need to convert FV and EV
    hpc2=STRUCT(
      MKPOLS(lar.V, lar.C[:FV]),
      MKPOLS(lar.V, lar.C[:EV])
    )

    VIEW(PROPERTIES(
        STRUCT(hpc1, T(2)(2.5),hpc2), 
        Properties("line_color" => WHITE,"line_width"=>2)
      ))
  end


end

# //////////////////////////////////////////////////////////////////////////////
function TestRandomLines()
  Random.seed!(0)
  hpc = STRUCT([RandomLine(2.0,3.0) for I in 1:6])
  lar = LAR(hpc)
  V,FVs,EVs = ARRANGE2D(lar.V,lar.C[:EV])
  VIEWCOMPLEX([Lar(V, Dict( :EV => EV)) for EV in EVs],explode=[1.2,1.2,1.2])
  VIEWCOMPLEX([Lar(V, Dict( :FV => FV)) for FV in FVs],explode=[1.2,1.2,1.2])
end


# ////////////////////////////////////////////////////////
function TestRandomBubbles()
  Random.seed!(0)
  hpc = STRUCT([RandomBubble() for I in 1:50])
  lar = LAR(hpc)
  V, EV  = lar.V, lar.C[:EV]
  V,FVs,EVs = ARRANGE2D(V,EV)
  VIEWCOMPLEX([Lar(V, Dict( :EV => EV)) for EV in EVs],explode=[1.2,1.2,1.2])
  VIEWCOMPLEX([Lar(V, Dict( :FV => FV)) for FV in FVs],explode=[1.2,1.2,1.2])
end

# //////////////////////////////////////////////////////////////////////////////
function TestTwoCubes()
  Random.seed!(0)
  cube = CUBOID([1,1,1])
  hpc = STRUCT(
    cube, 
    T(1,2,3)(.5,.5,0.75), 
    R(2,3)(pi/4), 
    R(1,3)(pi/4),
    cube)
  lar = LAR(hpc)
  
  arrangement = ARRANGE3D(lar)
  VIEWCOMPLEX(arrangement, explode=[1.2,1.2,2.0])
end

# ///////////////////////////////////////////////////////////
function TestSixCubes()
  Random.seed!(0)
  hpc = STRUCT([RandomCube(0.2,2.0) for I in 1:6])
  lar = LAR(hpc)
  arrangement = ARRANGE3D(lar)
  VIEWCOMPLEX(arrangement,explode=[1.4,1.4,1.4])
end

# ///////////////////////////////////////////////////////////
function TestCubeAndCylinders()
  Random.seed!(0)
  cyl = T(3)(-2)(CYLINDER([0.5,4])(8))
  hpc = STRUCT( 
    STRUCT(
      T(1,2,3)(-1,-1,-1),
      CUBE(2)), 
    cyl, R(2,3)(π/2), 
    cyl, R(1,3)(π/2), 
    cyl 
  )
  lar = LAR(hpc)
  arrangement = ARRANGE3D(lar)
  VIEWCOMPLEX(arrangement,explode=[1.2,1.2,1.8])
end


# /////////////////////////////////////////////////////////
function TestBool3D()

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
	
end

# //////////////////////////////////////////////////////////////////////////////
function TestLar()
  
  # TestToLAR()
  
  # TestRandomLines()
  # TestRandomBubbles()
  
  # TestTwoCubes()

  # TestSixCubes()
  TestCubeAndCylinders()

  # TestBool3D()

  println("TestLAR ok")
end

TestLar()

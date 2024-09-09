using Plasm

import Random

# ///////////////////////////////////////////////////////
function TestToLAR()

  # SPHERE: ::Hpc -> ::Lar
  obj = LAR(SPHERE(1.)([3,3]));
  V, FV, EV = obj.V, obj.C[:FV], obj.C[:EV];

  # SPHERE: ::Lar -> ::Hpc
  obj2 = MKPOLS(V,EV);
  VIEW(obj2)
  obj3 = MKPOLS(V,FV);
  VIEW(obj3)

  # SPHERE: ::Hpc -> ::Lar
  obj4 = LAR(obj3);
  V, FV, EV = obj4.V, obj4.C[:FV], obj4.C[:EV]

	# test2D
  if true
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
  if true
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
  if true
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

    VIEW(
      PROPERTIES(
        STRUCT(hpc1, T(2)(2.5),hpc2), 
        Properties("line_color" => WHITE,"line_width"=>2)
      ))
  end


end





# //////////////////////////////////////////////////////////////////////////////
function TestRandomLines()

  hpc = STRUCT([RandomLine(2.0,3.0) for I in 1:6])
  # VIEW(hpc)
  
  lar = LAR(hpc)
  V, EV  = lar.V, lar.C[:EV]
  V,FVs,EVs = arrange2d(V,EV)

  VIEW(COLORED(EXPLODECELLS(V, EVs, scale_factor=1.0)))
  VIEW(COLORED(EXPLODECELLS(V, EVs, scale_factor=1.2)))
  
  VIEW(COLORED(EXPLODECELLS(V, FVs, scale_factor=1.0)))
  VIEW(COLORED(EXPLODECELLS(V, FVs, scale_factor=1.2)))
  
end




# ////////////////////////////////////////////////////////
function TestRandomBubbles()

  hpc = STRUCT([RandomBubble() for I in 1:50])
  VIEW(hpc)
  
  lar = LAR(hpc)
  V, EV  = lar.V, lar.C[:EV]
  V,FVs,EVs = arrange2d(V,EV)
  
  VIEW(COLORED(EXPLODECELLS(V, EVs, scale_factor=1.0)))
  VIEW(COLORED(EXPLODECELLS(V, EVs, scale_factor=1.2)))
  VIEW(COLORED(EXPLODECELLS(V, FVs, scale_factor=1.0)))
  VIEW(COLORED(EXPLODECELLS(V, FVs, scale_factor=1.2)))
  
end


function TestRandomCubes()

  hpc = STRUCT([RandomCube(0.2,2.0) for I in 1:6])
  lar = LAR(hpc)
  V, copEV, copFE, copCF = arrange3d(lar)
  VIEWEXPLODED(TRIANGULATE3D(V, copEV, copFE, copCF)...)
  
end

# ///////////////////////////////////////////////////////////
function TestCubeAndCylinders()

  cyl = T(3)(-2)(CYLINDER([0.5,4])(4))

  hpc = STRUCT( 
    T(1,2,3)(-1,-1,-1)(CUBE(2)), 
    cyl, 
    R(2,3)(π/2), 
    cyl, 
    R(1,3)(π/2), 
    cyl 
  )
  VIEW(hpc)
  
  lar = LAR(hpc)
  V, copEV, copFE, copCF = arrange3d(lar)
  VIEWEXPLODED(TRIANGULATE3D(V, copEV, copFE, copCF)...)

end

# /////////////////////////////////////////////////////////
function TestBool3D()

	cube = CUBE(1);
	assembly = STRUCT(cube, 
		 R(1,3)(pi/4), T(1,2,3)(.5,.5,.5), cube, 
		 R(2,3)(pi/5), T(1,2,3)(.5,.5,.5), cube);
	
  lar=LAR(assembly)
  V_original=copy(lar.V)

	V, copEV, copFE, copCF = arrange3d(lar)
	VIEWEXPLODED(TRIANGULATE3D(V, copEV, copFE, copCF)...)

	EV = cop2lar(copEV)
	FE = cop2lar(copFE) 
	CF = cop2lar(copCF)
	boolmatrix = bool3d(assembly, V, EV, FE, CF)

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
	
	V,CVs,FVs,EVs = TRIANGULATE3D(V_original, copEV, copFE, copCF) 
	V,CVs,FVs,EVs = TRIANGULATE3D(subassembly3d(V_original, copEV, copFE, copCF, diff  )...)
	V,CVs,FVs,EVs = TRIANGULATE3D(subassembly3d(V_original, copEV, copFE, copCF, inters) ...)
	V,CVs,FVs,EVs = TRIANGULATE3D(subassembly3d(V_original, copEV, copFE, copCF, union )...)
	
	GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1))
	GL.VIEW(GL.GLExplode(V,EVs,1.,1.,1.,1,1))
  """
	
end

# //////////////////////////////////////////////////////////////////////////////
function TestLar()
  Random.seed!(0)
  
  # basic lar
  TestToLAR()
  
  TestRandomLines()
  TestRandomBubbles()
  TestRandomCubes()
  TestCubeAndCylinders()

  # BROKEN 
  # TestBool3D()

  println("TestLAR ok")
end

TestLar()

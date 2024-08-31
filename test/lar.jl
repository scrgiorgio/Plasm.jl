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
function TestTriangulation()
  (V, copEV, copFE, copCF) = (
    BYCOL([
    0.70710678118655 0.70710678118655 0.0; 
    2.829225697485796e-17 0.0 0.0; 
    -0.70710678118655 0.70710678118655 0.0; 
    0.0 1.4142135623731 0.0; 
    0.4142135623731001 1.0 0.0; 
    1.0683217716804808e-17 1.0 0.0; 
    0.70710678118655 0.70710678118655 1.0; 
    0.0 0.0 1.0; 0.0 1.4142135623731 1.0; 
    0.4142135623731001 1.0 1.0; 
    -0.70710678118655 0.70710678118655 1.0; 
    0.0 1.0 1.0; 1.0 0.0 0.0; 
    1.0 1.0 0.0; 1.0 0.0 1.0; 
    1.0 1.0 1.0]), 
    sparse(
      [1, 5, 10, 1, 2, 6, 9, 20, 2, 3, 16, 3, 4, 13, 4, 5, 7, 14, 21, 6, 7, 27, 8, 10, 11, 8, 9, 15, 18, 23, 12, 13, 17, 11, 12, 14, 19, 28, 15, 16, 17, 18, 19, 27, 20, 22, 24, 21, 22, 26, 23, 24, 25, 25, 26, 28], 
      [1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16], 
      Int8[-1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1],
      28, 16), 
    sparse(
      [2, 3, 10, 1, 6, 1, 7, 1, 4, 2, 5, 10, 1, 2, 13, 1, 2, 15, 3, 8, 16, 3, 6, 11, 13, 3, 5, 5, 8, 16, 4, 9, 4, 7, 4, 5, 14, 15, 6, 9, 6, 7, 7, 9, 8, 9, 13, 8, 9, 15, 10, 11, 10, 14, 10, 12, 11, 16, 11, 12, 12, 16, 12, 14, 13, 15, 14, 16 ], 
      [1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 10, 10, 11, 11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28], 
      Int8[-1, 1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, -1, 1, 1, -1], 
      16, 28), 
    sparse(
      [2, 4, 2, 3, 1, 3, 2, 4, 1, 3, 2, 4, 2, 4, 2, 3, 2, 4, 1, 2, 1, 2, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2], 
      [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16], 
      [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1], 
      4, 16))
  EV = AA(sort)([findnz(copEV[k,:])[1] for k=1:copEV.m]); # vertices per edge
  FE = AA(sort)([findnz(copFE[k,:])[1] for k=1:copFE.m]); # edges per face
  FV = [union(CAT([EV[e]  for e in f])) for f in FE]; # verts x face
  triangulated_faces = LAR2TRIANGLES(V, EV, FV, FE)
end   


# //////////////////////////////////////////////////////////////////////////////
function TestDrawAtoms()
  cube = CUBE(1)
  assembly = STRUCT(cube, R(1,2)(π/4), cube)
  lar=LAR(assembly)

  V,FV,EV = lar.V, lar.C[:FV], lar.C[:EV]	 

  copEV = cop_boundary_0(EV)
  copFE = cop_boundary_1(V, FV, EV) ## TODO: debug
  V, copEV, copFE, copCF = arrange3D(V, copEV, copFE )

  # generate and draw the atoms [and the b-rep of outer space]
  atoms,__CF = get_atoms(copEV,copFE,copCF)
  VIEWATOMS(V,copEV,copFE,copCF, atoms; view_outer=true)
end


# //////////////////////////////////////////////////////////////////////////////
function TestRandomLines()

  function RandomLine(size_min::Float64,size_max::Float64)
    size = size_min+rand()*(size_max-size_min)
    return STRUCT(
      T(1,2)(rand(2)...), 
      S([1,2])([size,size]), 
      R([1,2])(2*pi*rand()),
      Plasm.SQUARE(1)
    )
  end

  hpc = STRUCT([RandomLine(2.0,3.0) for I in 1:6])
  # VIEW(hpc)
  
  lar = LAR(hpc)
  V, EV  = lar.V, lar.C[:EV]
  V,FVs,EVs = arrange2D(V,EV)
  
  function ViewColored(V,cells, scale=1.2, line_width=3)
    exploded = explodecells(V, cells, sx=scale, sy=scale, sz=scale)
    v=[]
    for k in eachindex(exploded)
      c = Point4d(Plasm.COLORS[(k-1)%12+1] - (rand(Float64,4)*0.1))
      c[4] = 1.0
      push!(v,PROPERTIES(exploded[k], Properties(
        "line_color" => c, 
        "face_color" => c,
        "line_width" => line_width)))
    end
    VIEW(STRUCT(v))
  end


  ViewColored(V, EVs, 1.0)
  ViewColored(V, EVs, 1.2)
  
  ViewColored(V, FVs, 1.0)
  ViewColored(V, FVs, 1.2)
  
end



# //////////////////////////////////////////////////////////////////////////////
function TestRandomCubes()

  function RandomCube(size_min::Float64,size_max::Float64)
    size = size_min+rand()*(size_max-size_min)
    return STRUCT(
      T(1,2,3)(rand(3)...), 
      S([1,2,3])([size,size,size]), 
      R([1,2])(2*pi*rand()),
      R([2,3])(2*pi*rand()),
      R([1,3])(2*pi*rand()),
      Plasm.CUBE(1) 
    )
  end

  hpc = STRUCT([RandomCube(0.2,2.0) for I in 1:6])
  lar = LAR(hpc)
  
  V, EV, FV  = lar.V, lar.C[:EV], lar.C[:FV]

  copEV = cop_boundary_0(EV)
  copFE = cop_boundary_1(V, FV, EV)
  V, copEV, copFE, copCF = arrange3D(V, copEV, copFE)
  V,CVs,FVs,EVs = pols2tria(V, copEV, copFE, copCF);
  VIEWEXPLODED(V, CVs, FVs, EVs)
  
end


# ////////////////////////////////////////////////////////
function TestRandomBubbles()

  function RandomBubble()
    vs = rand()
    vt = rand(2)
    return STRUCT(
      T(1,2)(vt...),
      S([1,2])([0.25*vs, 0.25*vs]), 
      CIRCUMFERENCE(1)(rand(3:32))
    )
  end

  hpc = STRUCT([RandomBubble() for I in 1:50])
  VIEW(hpc)
  
  lar = LAR(hpc)
  V, EV  = lar.V, lar.C[:EV]
  V,FVs,EVs = arrange2D(V,EV)
  
  function ViewColored(V,cells, scale=1.2, line_width=3)
    exploded = explodecells(V, cells, sx=scale, sy=scale, sz=scale)
    v=[]
    for k in eachindex(exploded)
      c = Point4d(Plasm.COLORS[(k-1)%12+1] - (rand(Float64,4)*0.1))
      c[4] = 1.0
      push!(v,PROPERTIES(exploded[k], Properties(
        "line_color" => c, 
        "face_color" => c,
        "line_width" => line_width)))
    end
    VIEW(STRUCT(v))
  end

  ViewColored(V, EVs, 1.0)
  ViewColored(V, EVs, 1.2)
  ViewColored(V, FVs, 1.0)
  ViewColored(V, FVs, 1.2)
  
end


# /////////////////////////////////////////////////////////
function TestBool3D()

	cube = CUBE(1);
	assembly = STRUCT(cube, 
		 R(1,3)(pi/4), T(1,2,3)(.5,.5,.5), cube, 
		 R(2,3)(pi/5), T(1,2,3)(.5,.5,.5), cube);
	
  lar=LAR(assembly)
	V,FV,EV = new2old()
	
	copEV = cop_boundary_0(EV)
	copFE = cop_boundary_1(V, FV, EV)
  V_original=V
	V, copEV, copFE, copCF = arrange3D(V, copEV, copFE )

	boolmatrix = bool3d(assembly, V,copEV,copFE,copCF);
	V,CVs,FVs,EVs = pols2tria(V, copEV, copFE, copCF)
	VIEWEXPLODED(V,CVs,FVs,EVs)
	
	Matrix(boolmatrix)
	
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
	
	unione       = Matrix(copCF)' * Int.(AorBorC  ) # coord vector of Faces
	intersection = Matrix(copCF)' * Int.(AandBandC) # coord vector of Faces
	difference   = Matrix(copCF)' * Int.(AminBminC) # coord vector of Faces
	
	V,CVs,FVs,EVs = pols2tria(V_original, copEV, copFE, copCF) 
	V,CVs,FVs,EVs = pols2tria(V_original, copEV, copFE, copCF, difference) 
	V,CVs,FVs,EVs = pols2tria(V_original, copEV, copFE, copCF, intersection) 
	V,CVs,FVs,EVs = pols2tria(V_original, copEV, copFE, copCF, unione)
	
	GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1))
	GL.VIEW(GL.GLExplode(V,EVs,1.,1.,1.,1,1))
	
end

# //////////////////////////////////////////////////////////////////////////////
function TestLar()
  Random.seed!(0)
  
  # basic lar
  TestToLAR()
  
  # arrangements 2d
  TestRandomLines()
  TestRandomBubbles()

  # arrangmenets3d
  TestRandomCubes()

  # BROKEN
  # TestTriangulation()
  # TestDrawAtoms()
  # TestBool3D()

  println("TestLAR ok")
end

TestLar()

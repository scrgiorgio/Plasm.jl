using Plasm
using Test

# //////////////////////////////////////////////////////////////////////////////////////////
function TestComputeNormal()
	@test ComputeTriangleNormal([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]) == [0.0, 0.0, 1.0]
	@test ComputeTriangleNormal([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]) == [1.0, 0.0, 0.0]
	@test ComputeTriangleNormal([0.0, 0.0, 0.0], [1.0, 0.0, 1.0], [1.0, 0.0, 0.0]) == [0.0, 1.0, 0.0]
end

function TestGoodTet()
	@test GoodTetOrientation([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]) == true
	@test GoodTetOrientation([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]) == false
end

function TestBox()
	x1, y1, z1 = 1.0, 2.0, 3.0
	x2, y2, z2 = 4.0, 5.0, 6.0
	b = BoxNd([1.0, x1, y1, z1], [1.0, x2, y2, z2])
	@test dim(b) == 4
	@test b.p1 == [1.0, x1, y1, z1]
	@test b.p2 == [1.0, x2, y2, z2]
	@test valid(b) == true
	@test b == BoxNd(toList(b)...)
	@test size(b) == [0.0, x2 - x1, y2 - y1, z2 - z1]
	@test center(b) == [1.0, (x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2]
	x1, y1, z1 = -1, -2, -3
	x2, y2, z2 = 10, 20, 30
	addPoint(b, [1.0, x1, y1, z1])
	addPoint(b, [1.0, x2, y2, z2])
	@test b == BoxNd([1.0, x1, y1, z1], [1.0, x2, y2, z2])
end

function TestHpcMat()
  T = MatrixNd(4)
  v = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
  T[1, 1] += 0.0
  @test T == MatrixNd(v)
  @test dim(T) == 4
  @test T[1, 1] == 1.0
  @test isIdentity(T) == true
  @test toList(T) == v
  @test transpose(T) == T
  @test invert(T) == T
  @test embed(T, 5) == MatrixNd(5)
  @test T * T == T
  
  # transform point does not want homo coordinate
  p = [2.0, 3.0, 4.0]
  @test transformPoint(T, p) == p
  @test translate([10.0, 20.0, 30.0]) == MatrixNd([[1.0, 0.0, 0.0, 0.0], [10.0,  1.0, 0.0, 0.0], [20.0, 0.0,  1.0, 0.0], [30.0, 0.0, 0.0,  1.0]])
  @test scale([10.0, 20.0, 30.0])     == MatrixNd([[1.0, 0.0, 0.0, 0.0], [ 0.0, 10.0, 0.0, 0.0], [ 0.0, 0.0, 20.0, 0.0], [ 0.0, 0.0, 0.0, 30.0]])
  angle = π / 2
  @test rotate(1, 2, angle) == MatrixNd([[1.0, 0.0, 0.0], [0.0, cos(angle), -sin(angle)], [0.0, sin(angle), cos(angle)]])
end

function TestHpcMkPol()

  # 1D
  points=[[0.0],[1.0],[2.0],   [8.0],[9.0],[10.0]]
  hulls=[[1,2,3],[4,5,6]]
  obj=MkPol(points,hulls).childs[1]
  @test dim(obj)==1
  @test box(obj)==BoxNd([0.0],[10.0])
  obj=ToSimplicialForm(obj)
  @test length(obj.points)==4
  @test length(obj.hulls)==2
  @test box(obj)==BoxNd([0.0],[10.0])

  # 2D
  points = [[0.0, 0.0], [0.2, 0.2], [1.0, 0.0], [0.3, 0.3], [1.0, 1.0], [0.4, 0.4], [0.0, 1.0], [0.5, 0.5], [0.2, 0.8]]
  hulls = [collect(1:length(points))]
  obj = MkPol(points, hulls).childs[1]
  @test dim(obj) == 2
  @test box(obj) == BoxNd([0.0,0.0], [1.0,1.0])
  obj = ToSimplicialForm(obj)
  @test length(obj.points) == 4
  @test length(obj.hulls) == 2
  @test box(obj) == BoxNd([0.0,0.0], [1.0,1.0])

  # 3D
  points = [
	  [0.0, 0.0, 0.0],
	  [1.0, 0.0, 0.0],
	  [1.0, 1.0, 0.0],
	  [0.0, 1.0, 0.0],
	  [0.0, 0.0, 1.0],
	  [1.0, 0.0, 1.0],
	  [1.0, 1.0, 1.0],
	  [0.0, 1.0, 1.0],
	  [0.1, 0.1, 0.1],
	  [0.2, 0.2, 0.2],
	  [0.3, 0.3, 0.3]
  ]
  hulls = [collect(1:length(points))]
  obj = MkPol(points, hulls).childs[1]
  @test dim(obj) == 3
  @test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
  obj = ToSimplicialForm(obj)
  @test dim(obj) == 3
  @test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
  @test length(obj.points) == 8
  @test length(obj.hulls) == 6
end



function TestHpcInternal()

  points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.2, 0.2], [0.3, 0.3], [0.4, 0.4], [0.5, 0.5], [0.2, 0.8]]
  hulls = [collect(1:length(points))]
  obj = MkPol(points, hulls)
  @test dim(obj) == 2
  @test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

  # struct
  obj = Struct([obj])
  @test dim(obj) == 2
  @test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])
  
  # cube
  obj = Cube(3, 0.0, 1.0)
  @test dim(obj) == 3
  @test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
  @test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
  v = toList(obj)
  @test length(v) == 1
  (T, properties, obj) = v[1]
  obj = ToSimplicialForm(obj)
  @test dim(obj) == 3
  @test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
  @test length(obj.points) == 8
  @test length(obj.hulls) == 6
  
  # simplex
  obj = Simplex(3)
  @test dim(obj) == 3
  @test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
  v = toList(obj)
  @test length(v) == 1
  (T, properties, obj) = v[1]
  obj = ToSimplicialForm(obj)
  @test dim(obj) == 3
  @test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
  @test length(obj.points) == 4
  @test length(obj.hulls) == 1
  
  # quote
  obj = Quote([1.0, -2.0, 1.0])
  @test box(obj) == BoxNd([0.0], [4.0])
  (T, properties, obj) = toList(obj)[1]
  @test length(obj.points) == 4
  @test obj.hulls == [[1, 2], [3, 4]]
  
  # join
  obj = Join([Cube(2, 0.0, 1.0), Cube(2, 0.0, 1.0)])
  (T, properties, obj) = toList(obj)[1]
  obj = ToSimplicialForm(obj)
  @test length(obj.points) == 4
  @test length(obj.hulls) == 2
  
  # Transform
  obj = Transform(Cube(2, 0.0, 1.0),MatrixNd(3))
  (T, properties, obj) = toList(obj)[1]
  @test dim(T) == 3
  @test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

  # power
  obj = Power(Cube(2, 0.0, 1.0), Cube(1, 10.0, 20.0))
  l = toList(obj)
  @test length(l) == 1
  (T, properties, obj) = l[1]
  @test length(obj.points) == 8
  @test box(obj) == BoxNd([0.0, 0.0, 10.0], [1.0, 1.0, 20.0])

  # Scale
  obj = Scale(Cube(2, 0.0, 1.0),[1.0, 1.0])
  @test dim(obj) == 2
  @test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

  # Rotate
  obj = Rotate(Cube(2, 0.0, 1.0), 1,2, 0.0)
  @test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

  # Translate
  obj = Translate(Cube(2, 0.0, 1.0),[0.0, 0.0])
  @test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

  # MapFn
  function fn(p)
	  return [p[I]+0.0 for I in 1:length(p)]
  end
  obj=MapFn(Cube(2, 0.0, 1.0),fn)
  @test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

  # UkPol
  points,hulls=UkPol(Cube(2, 0.0, 1.0))
  @test length(points)==4
  @test length(hulls)==1 && length(hulls[1])==4

  # ToBoundaryForm
  obj=ToBoundaryForm(Struct([
	  Translate(Cube(2),[0.0,0.0]),
	  Translate(Cube(2),[1.0,0.0])
  ]))
  (T, properties, obj) = toList(obj)[1]
  @assert length(obj.hulls)==6
  @assert length(obj.points)==6

end

# ///////////////////////////////////////////////////////
function TestToLAR()

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

# ///////////////////////////////////////////////////////
function MyMain()
  TestComputeNormal()
  TestGoodTet()
  TestBox()
  TestHpcMat()
  TestHpcMkPol()
  TestHpcInternal()
  TestToLAR()
  println("TestHpc ok")
end


MyMain()



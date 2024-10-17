using Plasm
using Test

# ////////////////////////////////////////////////////////////////
function TestViewer()

  begin
    viewer=Viewer()
    render_cuboid(viewer, Box3d(Point3d(0, 0, 0), Point3d(1, 1, 1)))
    render_axis(viewer, Point3d(0, 0, 0), Point3d(+1.1, +1.1, +1.1))
    run_viewer(viewer)
  end

  begin
    viewer=Viewer()
    render_text(viewer, "hello")
    run_viewer(viewer)
  end

  begin

    viewer=Viewer()
    
    render_points(viewer, Vector{Float32}([0.0,0.0,0.0, 1.0,1.0,0.0 ]); point_size=10, point_color=YELLOW)
    render_points(viewer, Vector{Float32}([0.0,0.0,0.2, 1.0,1.0,0.2 ]),colors=Vector{Float32}([1.0,0.0,0.0,1.0,  0.0,1.0,0.0,1.0]); point_size=10)
    
    render_lines(viewer, Vector{Float32}([0.0,0.0,0.4, 1.0,1.0,0.4 ]); line_width=10, line_color=YELLOW)
    render_lines(viewer, Vector{Float32}([0.0,0.0,0.6, 1.0,1.0,0.6 ]); colors=Vector{Float32}([1.0,0.0,0.0,1.0,  0.0,1.0,0.0,1.0]),line_width=10)
    
    render_triangles(viewer, Vector{Float32}([0.0,0.0,0.8, 1.0,0.0,0.8, 1.0,1.0,0.8 ]); face_color=YELLOW)
    render_triangles(viewer, Vector{Float32}([0.0,0.0,1.0, 1.0,0.0,1.0, 1.0,1.0,1.0 ]); colors=Vector{Float32}([1.0,0.0,0.0,1.0,  0.0,1.0,0.0,1.0,  0.0,0.0,1.0,1.0]))
    
    run_viewer(viewer)

  end


end

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
  @test translate([10.0, 20.0, 30.0]) == MatrixNd([[1.0, 0.0, 0.0, 0.0], [10.0, 1.0, 0.0, 0.0], [20.0, 0.0, 1.0, 0.0], [30.0, 0.0, 0.0, 1.0]])
  @test scale([10.0, 20.0, 30.0]) == MatrixNd([[1.0, 0.0, 0.0, 0.0], [0.0, 10.0, 0.0, 0.0], [0.0, 0.0, 20.0, 0.0], [0.0, 0.0, 0.0, 30.0]])
  angle = π / 2
  @test rotate(1, 2, angle) == MatrixNd([[1.0, 0.0, 0.0], [0.0, cos(angle), -sin(angle)], [0.0, sin(angle), cos(angle)]])
end

function TestHpcMkPol()

  # 1D
  points = [[0.0], [1.0], [2.0], [8.0], [9.0], [10.0]]
  hulls = [[1, 2, 3], [4, 5, 6]]
  obj = MkPol(points, hulls).childs[1]
  @test dim(obj) == 1
  @test box(obj) == BoxNd([0.0], [10.0])
  obj = ToSimplicialForm(obj)
  @test length(obj.points) == 4
  @test length(obj.hulls) == 2
  @test box(obj) == BoxNd([0.0], [10.0])

  # 2D
  points = [[0.0, 0.0], [0.2, 0.2], [1.0, 0.0], [0.3, 0.3], [1.0, 1.0], [0.4, 0.4], [0.0, 1.0], [0.5, 0.5], [0.2, 0.8]]
  hulls = [collect(1:length(points))]
  obj = MkPol(points, hulls).childs[1]
  @test dim(obj) == 2
  @test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])
  obj = ToSimplicialForm(obj)
  @test length(obj.points) == 4
  @test length(obj.hulls) == 2
  @test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

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


# /////////////////////////////////////////////////////////////////////////
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

  # HpcSimplex
  obj = HpcSimplex(3)
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
  obj = Transform(Cube(2, 0.0, 1.0), MatrixNd(3))
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
  obj = Scale(Cube(2, 0.0, 1.0), [1.0, 1.0])
  @test dim(obj) == 2
  @test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

  # Rotate
  obj = Rotate(Cube(2, 0.0, 1.0), 1, 2, 0.0)
  @test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

  # Translate
  obj = Translate(Cube(2, 0.0, 1.0), [0.0, 0.0])
  @test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

  # MapFn
  function fn(p)
    return [p[I] + 0.0 for I in 1:length(p)]
  end
  obj = MapFn(Cube(2, 0.0, 1.0), fn)
  @test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

  # UkPol
  points, hulls = UkPol(Cube(2, 0.0, 1.0))
  @test length(points) == 4
  @test length(hulls) == 1 && length(hulls[1]) == 4

  # ToBoundaryForm
  obj = ToBoundaryForm(Struct([
    Translate(Cube(2), [0.0, 0.0]),
    Translate(Cube(2), [1.0, 0.0])
  ]))
  (T, properties, obj) = toList(obj)[1]
  @assert length(obj.hulls) == 6
  @assert length(obj.points) == 6

end


# ///////////////////////////////////////////////////////////////////////
function Temple()

  function Column(r::Float64, h::Float64)
    basis = CUBOID([2 * r * 1.2, 2 * r * 1.2, h / 12.0])
    trunk = CYLINDER([r, (10.0 / 12.0) * h])(12)
    capital = basis
    beam = S(1)(3)(capital)
    return TOP([TOP([TOP([basis, trunk]), capital]), beam])
  end

  function Gable(radius::Float64, h::Float64, n::Int64)
    lastX = n * 3 * (2 * radius * 1.2)
    triangle = MKPOL(
      [[0, 0], [lastX, 0], [lastX / 2.0, h / 2]],
      [[1, 2, 3]],
      [[1]])
    return R([2, 3])(PI / 2)(POWER([triangle, QUOTE([radius * 1.2])]))
  end

  col = Column(1.0, 12.0)

  function ColRow(N::Int64)
    v = [col for I in 1:N]
    return INSR(RIGHT)(v)
  end

  ColRowAndGable = TOP([ColRow(4), Gable(1.0, 12.0, 4)])

  Temple = STRUCT(CAT([
    [ColRowAndGable, T(2)(6.0)],
    DOUBLE_DIESIS(4)([ColRow(4), T(2)(6.0)]),
    [ColRowAndGable]
  ]))

  Ground = EMBED(1)(BOX([1, 2])(Temple))

  Xsizes = QUOTE(DOUBLE_DIESIS(12)([0.6, -1.2]))
  Ysizes = QUOTE(DOUBLE_DIESIS(5)([-1, 5]))
  Zsizes = QUOTE([-13, 0.6])

  SecondaryBeams = POWER([POWER([Xsizes, Ysizes]), Zsizes])
  model = STRUCT([Temple, Ground, T(1, 2, 3)(-1.65, 0.6, 0), SecondaryBeams,])
  VIEW(model)

  skel1model = SKELETON(1)(model)
  VIEW(skel1model)
  skel1model = SKELETON(2)(model)
  VIEW(skel1model)

  # SKELETON(4)(model) -> error
end

import Base.*, Base.splat
*(a::Hpc, b::Hpc) = Power(a, b)

function Manhattan2D()
  global verts = [[0.0, 0], [3, 0], [5, 0], [7, 0], [8, 0], [9.5, 1], [10, 1.5], [0, 3], [3, 3],
    [5, 3], [7, 3], [8, 3], [9.5, 3], [0, 4], [3, 4], [5, 4], [9.5, 4], [12, 4], [9.5, 5], [10, 5],
    [12, 5], [0, 6], [3, 6], [5, 6], [0, 7], [3, 7], [5, 7], [9.5, 7], [12, 7], [9.5, 8],
    [12, 8], [0, 9], [3, 9], [5, 9], [8, 9], [9, 9], [12, 9], [0, 10], [3, 10], [5, 10],
    [8, 10], [9, 10], [9.5, 10], [10, 10], [12, 10], [6, 11], [7, 11], [0, 12], [3, 12], [9, 12],
    [9.5, 12], [0, 13], [3, 13], [6, 13], [7, 13], [9, 13], [9.5, 13], [0, 14], [3, 14], [5, 14],
    [8, 14], [9, 14], [9.5, 14], [10, 14], [12, 14], [0, 15], [3, 15], [5, 15], [8, 15], [0, 16],
    [6, 16], [7, 16], [9, 17], [9.5, 17], [10, 17], [12, 17], [6, 18], [7, 18], [9, 18], [9.5, 18],
    [10, 18], [12, 18], [2, 19], [3, 19], [5, 19], [8, 19], [9, 19], [9.5, 19], [10, 19], [12, 19],
    [5, 20], [12, 20], [7, 22], [10, 22], [9, 6], [12, 6], [9, 15], [9.5, 15], [10, 15], [12, 15]]
  global cells = [[1, 2, 9, 8], [3, 4, 11, 10], [5, 6, 13, 12], [14, 15, 23, 22], [16, 17, 19, 24],
    [7, 18, 21, 20], [25, 26, 33, 32], [27, 95, 28, 35, 34], [95, 96, 29, 28], [30, 31, 37, 36],
    [38, 39, 49, 48], [40, 41, 47, 46], [41, 61, 55, 47], [55, 61, 60, 54], [54, 60, 40, 46],
    [42, 43, 51, 50], [44, 45, 65, 64], [52, 53, 59, 58], [56, 57, 63, 62], [66, 67, 84, 83, 70],
    [68, 69, 72, 71], [69, 86, 78, 72], [78, 86, 85, 77], [71, 77, 85, 68], [97, 98, 74, 73],
    [99, 100, 76, 75], [79, 80, 88, 87], [81, 82, 90, 89], [91, 92, 94, 93]]
  model = MKPOL(verts, cells)
  VIEW(model, title="Manhattan2D")
end

function Manhattan3D()
  ManhattanH = [1, 3, 1, 11, 1, 2, 1, 1, 1, 8, 15, 1, 1, 1, 1, 8, 1, 15, 8, 1, 2, 2, 2, 2, 5, 9, 1, 1, 1] .* 3
  storeys = CONS(AA(DIESIS)(ManhattanH))(0.5)
  pols1D = AA(QUOTE)(storeys)
  pols2D = [MKPOL(verts, [cell]) for cell in cells]
  pols3D = AA(splat(*))(TRANS([pols2D, pols1D]))
  VIEW(STRUCT(pols3D), title="Manhattan3D")
end


function TestProperties()

  # example: how to set the overall background color
  VIEW(CUBE(1),
  properties=Properties(
      "background_color" => DEFAULT_BACKGROUND_COLOR,
      "title" => "cube with background color"
    ))

  # line_with property
  VIEW(
    PROPERTIES(
      CUBE(1),
      Properties("line_width" => 10)),
    title="line_width")

  # line_color property
  VIEW(
    PROPERTIES(
      CUBE(1),
      Properties("line_color" => Point4d(1.0, 0.0, 0.0, 1.0))),
    title="line_color"
  )

  # face_color property
  VIEW(
    PROPERTIES(
      CUBE(1),
      Properties("face_color" => Point4d(0.0, 1.0, 0.0, 1.0))),
    title="face_color"
  )

  # COLOR
  VIEW(COLOR([1.0, 1.0, 0.0, 0.0])(CUBE(1)))

  # example of mixing and matching different properties
  cube_a = STRUCT(T(1)(0.0), CUBE(1))
  a = PROPERTIES(cube_a,
    Properties(
      "face_color" => Point4d(0.0, 1.0, 0.0, 1.0),
      "line_color" => Point4d(1.0, 1.0, 0.0, 1.0),
      "line_width" => DEFAULT_LINE_WIDTH
    )
  )
  cube_b = STRUCT(T(1)(1.0), CUBE(1))
  b = PROPERTIES(cube_b,
    Properties(
      "face_color" => Point4d(1.0, 0.0, 0.0, 1.0),
      "line_color" => Point4d(0.0, 1.0, 1.0, 1.0),
      "line_width" => DEFAULT_LINE_WIDTH
    )
  )
  VIEW(
    STRUCT(a, b),
    properties=Properties(
      "title" => "2 colored cubes",
      "background_color" => DEFAULT_BACKGROUND_COLOR
    )
  )

  # //////////////////////////////////////////////////////
  # 2d frame/hull/point example
  begin

    points = [[rand() + 1.0, rand() + 1.0] for I in 1:100]

    obj = STRUCT(

      # show points
      PROPERTIES(
        MKPOINTS(points),
        Properties(
          "point_color" => YELLOW,
          "point_size" => 3
        )
      ),

      # show hull (note face color is transparent, so no fill)
      PROPERTIES(
        MKPOL(
          points,
          [[it for it in 1:length(points)]]
        ),
        Properties(
          "face_color" => TRANSPARENT,
          "line_color" => GREEN,
          "line_width" => DEFAULT_LINE_WIDTH
        )
      ),
      # show frame
      FRAME2([0.0, 0.0], [1.0, 1.0]),
    )

    VIEW(obj, properties=Properties("show_axis" => false))
  end

  # //////////////////////////////////////////////////////
  # 3d frame/hull/point example
  begin

    points = [[rand() + 1.0, rand() + 1.0, rand() + 1.0] for I in 1:100]

    obj = STRUCT(

      # show points
      PROPERTIES(
        MKPOINTS(points),
        Properties(
          "point_color" => YELLOW,
          "point_size" => 3
        )
      ),
      # show hull
      PROPERTIES(
        MKPOL(
          points,
          [[it for it in 1:length(points)]]
        ),
        Properties(
          "face_color" => GRAY,
          "line_color" => GREEN,
          "line_width" => DEFAULT_LINE_WIDTH
        )
      ),
      # show frame
      FRAME3(Point3d(0, 0, 0), Point3d(1, 1, 1)),
    )

    VIEW(obj, properties=Properties("show_axis" => false))

  end
end


# ////////////////////////////////////////////////
function TestProperties2()

  # ////////////////////////////////////////////////////////
  # geometry (per-batch) properties 
  # ///////////////////////////////////////////////////////
  obj = PROPERTIES(
    CUBOID([1, 1, 1]),
    Properties(
      "point_size" => DEFAULT_POINT_SIZE,
      "point_color" => DEFAULT_POINT_COLOR,
      "line_width" => DEFAULT_LINE_WIDTH,
      "line_color" => DEFAULT_LINE_COLOR,
      "face_color" => DEFAULT_FACE_COLOR,

      # used by certain algorithms to set/retrieve node type
      "node_type" => "whatever-you-need",
    ))

  # ////////////////////////////////////////////////////////
  # viewer (global) properties
  # ////////////////////////////////////////////////////////
  VIEW(obj, properties=Properties(
    "title" => "My viewer example",
    "background_color" => DEFAULT_BACKGROUND_COLOR,
  ))

end


# /////////////////////////////////////////////////////////////
function TestSphere()
  VIEW(SPHERE(1.0)([16, 16]), title="TestSphere")
end

# /////////////////////////////////////////////////////////////
function TestTorus()
  VIEW(TORUS([1.0, 2.0])([20, 20]), title="TestTorus")
end

# /////////////////////////////////////////////////////////////
function TestBezier()
  VIEW(MAP(BEZIER(S1)([[-0, 0], [1, 0], [1, 1], [2, 1], [3, 1]]))(INTERVALS(1.0)(32)), title="TestBezier-1")
  C0 = BEZIER(S1)([[0, 0, 0], [10, 0, 0]])
  C1 = BEZIER(S1)([[0, 2, 0], [8, 3, 0], [9, 2, 0]])
  C2 = BEZIER(S1)([[0, 4, 1], [7, 5, -1], [8, 5, 1], [12, 4, 0]])
  C3 = BEZIER(S1)([[0, 6, 0], [9, 6, 3], [10, 6, -1]])
  VIEW(MAP(BEZIER(S2)([C0, C1, C2, C3]))(Power(INTERVALS(1.0)(10), INTERVALS(1.0)(10))), title="TestBezier-2")
end

function TestCoonsPatch()
  Su0 = BEZIER(S1)([[0, 0, 0], [10, 0, 0]])
  Su1 = BEZIER(S1)([[0, 10, 0], [2.5, 10, 3], [5, 10, -3], [7.5, 10, 3], [10, 10, 0]])
  Sv0 = BEZIER(S2)([[0, 0, 0], [0, 0, 3], [0, 10, 3], [0, 10, 0]])
  Sv1 = BEZIER(S2)([[10, 0, 0], [10, 5, 3], [10, 10, 0]])
  VIEW(MAP(COONSPATCH([Su0, Su1, Sv0, Sv1]))(Power(INTERVALS(1.0)(10), INTERVALS(1.0)(10))), title="TestCoonsPatch")
end

# /////////////////////////////////////////////////////////////
function TestRuledSurface()
  alpha = point -> [point[1], point[1], 0]
  beta = point -> [-1, +1, point[1]]
  domain = T([1, 2])([-1, -1])(Power(INTERVALS(2.0)(10), INTERVALS(2.0)(10)))
  VIEW(MAP(RULEDSURFACE([alpha, beta]))(domain), title="TestRuledSurface")
end

function TestProfileProdSurface()
  alpha = BEZIER(S1)([[0.1, 0, 0], [2, 0, 0], [0, 0, 4], [1, 0, 5]])
  beta = BEZIER(S2)([[0, 0, 0], [3, -0.5, 0], [3, 3.5, 0], [0, 3, 0]])
  domain = Power(INTERVALS(1.0)(20), INTERVALS(1.0)(20))
  VIEW(Struct([MAP(alpha)(domain), MAP(beta)(domain), MAP(PROFILEPRODSURFACE([alpha, beta]))(domain)]), title="TestProfileProdSurface")
end

function TestRotationalSurface()
  profile = BEZIER(S1)([[0, 0, 0], [2, 0, 1], [3, 0, 4]])
  domain = Power(INTERVALS(1.0)(10), INTERVALS(2 * PI)(30))
  VIEW(MAP(ROTATIONALSURFACE(profile))(domain), title="TestRotationalSurface")
end

DISPLAYNUBSPLINE
function TestConicalSurface()
  domain = Power(INTERVALS(1.0)(20), INTERVALS(1.0)(6))
  beta = BEZIER(S1)([[1, 1, 0], [-1, 1, 0], [1, -1, 0], [-1, -1, 0]])
  VIEW(MAP(CONICALSURFACE([[0, 0, 1], beta]))(domain), title="TestConicalSurface")
end

function TestCubicHermite()
  domain = INTERVALS(1.0)(20)
  VIEW(Struct([
      MAP(CUBICHERMITE(S1)([[1, 0], [1, 1], [-1, 1], [1, 0]]))(domain),
      MAP(CUBICHERMITE(S1)([[1, 0], [1, 1], [-2, 2], [2, 0]]))(domain),
      MAP(CUBICHERMITE(S1)([[1, 0], [1, 1], [-4, 4], [4, 0]]))(domain),
      MAP(CUBICHERMITE(S1)([[1, 0], [1, 1], [-10, 10], [10, 0]]))(domain)]), title="TestCubicHermite-1")
  c1 = CUBICHERMITE(S1)([[1, 0, 0], [0, 1, 0], [0, 3, 0], [-3, 0, 0]])
  c2 = CUBICHERMITE(S1)([[0.5, 0, 0], [0, 0.5, 0], [0, 1, 0], [-1, 0, 0]])
  sur3 = CUBICHERMITE(S2)([c1, c2, [1, 1, 1], [-1, -1, -1]])
  domain = Power(INTERVALS(1.0)(14), INTERVALS(1.0)(14))
  VIEW(MAP(sur3)(domain), title="TestCubicHermite-2")
end

function TestPermutahedron()
  VIEW(Struct([PERMUTAHEDRON(2), (PERMUTAHEDRON(2))]), title="TestPermutahedron-1")
  VIEW(Struct([PERMUTAHEDRON(3), (PERMUTAHEDRON(3))]), title="TestPermutahedron-2")
end

function TestSchegel3d()
  VIEW(SCHLEGEL3D(0.2)((T([1, 2, 3, 4])([-1.0 / 3.0, -1.0 / 3.0, -1, +1])(HPCSIMPLEX(4)))), title="TestSchegel3d-1")
  VIEW(SCHLEGEL3D(0.2)((T([1, 2, 3, 4])([-1, -1, -1, 1])(CUBOID([2, 2, 2, 2])))), title="TestSchegel3d-2")
  VIEW(SCHLEGEL3D(0.2)((T([1, 2, 3, 4])([-1.0 / 3.0, -1.0 / 3.0, -1, +1])(Power(HPCSIMPLEX(2), HPCSIMPLEX(2))))), title="TestSchegel3d-3")
end

function TestCubicSpline()
  domain = INTERVALS(1.0)(20)
  points = [
    [-3.0, 6.0],
    [-4.0, 2.0],
    [-3.0, -1.0],
    [-1.0, 1.0],
    [1.5, 1.5],
    [3.0, 4.0],
    [5.0, 5.0],
    [7.0, 2.0],
    [6.0, -2.0],
    [2.0, -3.0]
  ]
  VIEW(SPLINE(CUBICCARDINAL(domain))(points), title="TestCubicSpline-1")
  VIEW(SPLINE(CUBICUBSPLINE(domain))(points), title="TestCubicSpline-2")
end

function TestBilinarSurface()
  controlpoints = [
    [[0.0, 0.0, 0.0], [2.0, -4.0, 2.0]],
    [[0.0, 3.0, 1.0], [4.0, 0.0, 0.0]]
  ]
  domain = Power(INTERVALS(1.0)(10), INTERVALS(1.0)(10))
  mapping = BILINEARSURFACE(controlpoints)
  VIEW(MAP(mapping)(domain), title="TestBilinarSurface")
end

function TestBiquadraticSurface()
  controlpoints = [[[0, 0, 0], [2, 0, 1], [3, 1, 1]], [[1, 3, -1], [3, 2, 0], [4, 2, 0]], [[0, 9, 0], [2, 5, 1], [3, 3, 2]]]
  domain = Power(INTERVALS(1.0)(10), INTERVALS(1.0)(10))
  mapping = BIQUADRATICSURFACE(controlpoints)
  VIEW(MAP(mapping)(domain), title="TestBiquadraticSurface")
end



function TestHermiteSurface()
  controlpoints = ToFloat64([[[0, 0, 0], [2, 0, 1], [3, 1, 1], [4, 1, 1]], [[1, 3, -1], [3, 2, 0], [4, 2, 0], [4, 2, 0]], [[0, 4, 0], [2, 4, 1], [3, 3, 2], [5, 3, 2]], [[0, 6, 0], [2, 5, 1], [3, 4, 1], [4, 4, 0]]])
  domain = Power(INTERVALS(1.0)(10), INTERVALS(1.0)(10))
  mapping = HERMITESURFACE(controlpoints)
  VIEW(MAP(mapping)(domain), title="TestHermiteSurface")
end

function TestBezierSurface()
  controlpoints = [[[0, 0, 0], [0, 3, 4], [0, 6, 3], [0, 10, 0]], [[3, 0, 2], [2, 2.5, 5], [3, 6, 5], [4, 8, 2]], [[6, 0, 2], [8, 3, 5], [7, 6, 4.5], [6, 10, 2.5]], [[10, 0, 0], [11, 3, 4], [11, 6, 3], [10, 9, 0]]]
  domain = Power(INTERVALS(1.0)(10), INTERVALS(1.0)(10))
  mapping = BEZIERSURFACE(controlpoints)
  VIEW(MAP(mapping)(domain), title="TestBezierSurface")
end

function TestBezierManifold()
  grid1D = INTERVALS(1.0)(5)
  domain3D = Power(Power(grid1D, grid1D), grid1D)
  degrees = [2, 2, 2]
  Xtensor = [[[0, 1, 2], [-1, 0, 1], [0, 1, 2]], [[0, 1, 2], [-1, 0, 1], [0, 1, 2]], [[0, 1, 2], [-1, 0, 1], [0, 1, 2]]]
  Ytensor = [[[0, 0, 0.8], [1, 1, 1], [2, 3, 2]], [[0, 0, 0.8], [1, 1, 1], [2, 3, 2]], [[0, 0, 0.8], [1, 1, 1], [2, 3, 2]]]
  Ztensor = [[[0, 0, 0], [0, 0, 0], [0, 0, 0]], [[1, 1, 1], [1, 1, 1], [1, 1, 1]], [[2, 2, 1], [2, 2, 1], [2, 2, 1]]]
  mapping = BEZIERMANIFOLD(degrees)([Xtensor, Ytensor, Ztensor])
  VIEW(MAP(mapping)(domain3D), title="TestBezierManifold")
end

function TestOffset()
  verts = ToFloat64([[0, 0, 0], [3, 0, 0], [3, 2, 0], [0, 2, 0], [0, 0, 1.5], [3, 0, 1.5], [3, 2, 1.5], [0, 2, 1.5], [0, 1, 2.2], [3, 1, 2.2]])
  cells = [[1, 2], [2, 3], [3, 4], [4, 1], [5, 6], [6, 7], [7, 8], [8, 5], [1, 5], [2, 6], [3, 7], [4, 8], [5, 9], [8, 9], [6, 10], [7, 10], [9, 10]]
  pols = [[1]]
  House = MKPOL(verts, cells, pols)
  VIEW(Struct([OFFSET([0.1, 0.2, 0.1])(House), T(1)(1.2 * SIZE(1)(House))(House)]), title="TestOffset")
end

function TestThinSolid()
  Su0 = COMP([BEZIERCURVE([[0, 0, 0], [10, 0, 0]]), CONS([S1])])
  Su1 = COMP([BEZIERCURVE([[0, 10, 0], [2.5, 10, 3], [5, 10, -3], [7.5, 10, 3], [10, 10, 0]]), CONS([S1])])
  S0v = COMP([BEZIERCURVE([[0, 0, 0], [0, 0, 3], [0, 10, 3], [0, 10, 0]]), CONS([S2])])
  S1v = COMP([BEZIERCURVE([[10, 0, 0], [10, 5, 3], [10, 10, 0]]), CONS([S2])])
  surface = COONSPATCH([Su0, Su1, S0v, S1v])
  VIEW(MAP(surface)(Power(INTERVALS(1.0)(10), INTERVALS(1.0)(10))), title="TestThinSolid-1")
  solidMapping = THINSOLID(surface)
  Domain3D = Power(Power(INTERVALS(1.0)(5), INTERVALS(1.0)(5)), INTERVALS(0.5)(5))
  VIEW(MAP(solidMapping)(Domain3D), title="TestThinSolid-2")
end

function TestEllipse()
  VIEW(ELLIPSE([1.0, 2.0])(8), title="TestEllipse")
end

function TestBezierStripe()
  vertices = ToFloat64([[0, 0], [1.5, 0], [-1, 2], [2, 2], [2, 0]])
  VIEW(Struct([POLYLINE(vertices), Power(BEZIERSTRIPE([vertices, 0.25, 22]), QUOTE([0.9]))]), title="TestBezierStripe")
end

function TestDisplayNubSpline()
  degree = 3
  ControlPoints = ToFloat64([[0, 0], [-1, 2], [1, 4], [2, 3], [1, 1], [1, 2], [2.5, 1], [2.5, 3], [4, 4], [5, 0]]) # 10 points
  knots = [0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7] # 14 knots
  VIEW(DISPLAYNUBSPLINE(degree, knots, ControlPoints), title="TestDisplayNubSpline")
end

function TestDisplayNurbsSpline()
  knots = [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]
  _p = sqrt(2) / 2.0
  controlpoints = ToFloat64([[-1, 0, 1], [-_p, _p, _p], [0, 1, 1], [_p, _p, _p], [1, 0, 1], [_p, -_p, _p], [0, -1, 1], [-_p, -_p, _p], [-1, 0, 1]])
  VIEW(DISPLAYNURBSPLINE([2, knots, controlpoints]), title="TestDisplayNurbsSpline")
end

function TestMinkowski()
  p = MKPOL(
    [[0.0, 0.0]],
    [[1]],
    [[1]])
  B = MINKOWSKI([[-1.0 / 2.0, -1 * sqrt(3.0 / 2.0)], [-1.0 / 2.0, sqrt(3.0 / 2.0)], [1, 0]])(p)
  vertices = [[0, 0], [1, 0], [1, 0.5], [0.5, 0.5], [0.5, 1], [0, 1]]
  pol1D = MKPOL(vertices, [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]], [[1], [2], [3], [4], [5], [6]])
  pol2D = MKPOL(vertices, [[1, 2, 3, 4], [4, 5, 6, 1]], [[1, 2]])
  Min0 = STRUCT([T([1, 2])(v)(S([1, 2])([0.1, 0.1])(B)) for v in vertices])
  Min1 = MINKOWSKI([[0.1 * -1.0 / 2.0, 0.1 * -1 * sqrt(3.0 / 2.0)], [0.1 * -1.0 / 2.0, 0.1 * sqrt(3.0 / 2.0)], [0.1 * 1, 0.1 * 0]])(pol1D)
  Min2 = MINKOWSKI([[0.1 * -1.0 / 2.0, 0.1 * -1 * sqrt(3.0 / 2.0)], [0.1 * -1.0 / 2.0, 0.1 * sqrt(3.0 / 2.0)], [0.1 * 1, 0.1 * 0]])(pol2D)
  VIEW(Struct([Min0, Min1, Min2]), title="TestMinkowski")
end

function TestThinsolid()
  Su0 = COMP([BEZIERCURVE([[0, 0, 0], [10, 0, 0]]), CONS([S1])])
  Su1 = COMP([BEZIERCURVE([[0, 10, 0], [2.5, 10, 3], [5, 10, -3], [7.5, 10, 3], [10, 10, 0]]), CONS([S1])])
  S0v = COMP([BEZIERCURVE([[0, 0, 0], [0, 0, 3], [0, 10, 3], [0, 10, 0]]), CONS([S2])])
  S1v = COMP([BEZIERCURVE([[10, 0, 0], [10, 5, 3], [10, 10, 0]]), CONS([S2])])
  surface = COONSPATCH([Su0, Su1, S0v, S1v])
  VIEW(MAP(surface)(Power(INTERVALS(1.0)(10), INTERVALS(1.0)(10))), title="TestThinsolid-1")
  solidMapping = THINSOLID(surface)
  Domain3D = Power(Power(INTERVALS(1.0)(5), INTERVALS(1.0)(5)), INTERVALS(0.5)(5))
  VIEW(MAP(solidMapping)(Domain3D), title="TestThinsolid-2")
end

function TestSOLIDHELICOID(; nturns=3, R=1.0, r=0.0, shape=[36 * nturns, 8], pitch=2, thickness=0.1)
  totalangle = nturns * 2 * pi
  grid2D = Power(INTERVALS(36 * nturns)(36 * nturns), INTERVALS(4)(8))
  Domain2D = T([2])([r])(S([1, 2])([totalangle / shape[1], R - r])(grid2D))
  surface = p -> let (u, v) = p
    [v * cos(u); v * sin(u); u * (pitch / (2 * pi))]
  end
  VIEW(MAP(surface)(Domain2D))
  solidMapping = THINSOLID(surface)
  Domain3D = Power(Domain2D, INTERVALS(thickness)(1))
  VIEW(MAP(solidMapping)(Domain3D))
end;




function TestPolar()
  VIEW(POLAR(CUBOID([1, 1, 1])), title="TestPolar")
end

function TestSolidify()
  VIEW(SOLIDIFY(STRUCT(AA(POLYLINE)([
      [[0, 0], [4, 2], [2.5, 3], [4, 5], [2, 5], [0, 3], [-3, 3], [0, 0]],
      [[0, 3], [0, 1], [2, 2], [2, 4], [0, 3]],
      [[2, 2], [1, 3], [1, 2], [2, 2]]]))), title="TestSolidify")
end

function TestDiff()
  mypol1 = T([1, 2])([-5, -5])(CUBOID([10, 10]))
  mypol2 = S([1, 2])([0.9, 0.9])(mypol1)
  mypol3 = DIFF([mypol1, mypol2])
  VIEW(STRUCT([
      EX([0, 10])(mypol3), T(1)(12),
      LEX([0, 10])(mypol3), T(1)(25),
      S(3)(3)(SEX([0, PI])(16)(mypol3))
    ]), "TestDiff")
end

function TestCube()
  VIEW(Cube(3), title="TestCube")
end

function TestMapSphere()
  N, M = 16, 16
  domain = Power(
    Translate(Quote([pi / N for I in 1:N]), [-pi / 2]),
    Quote([2 * pi / M for I in 1:M]))
  fn = p -> [
    cos(p[1]) * sin(p[2]),
    cos(p[1]) * cos(p[2]),
    sin(p[1])
  ]
  obj = MAP(fn)(domain)
  VIEW(obj, title="TestMapSphere")
end

function TestFenvsMkPol()
  out = MkPol(
    [[0.0], [1.0], [2.0], [3.0], [4.0], [5.0]],
    [[6, 4], [1, 2]])
  VIEW(out, title="TestFenvsMkPol")
end

function TestArguments()

  STRUCT([CUBE(1), CUBE(1)])
  STRUCT(CUBE(1), CUBE(1))

  T(1, 2)(0.0, 0.0)
  T([1, 2])([0.0, 0.0])

  S(1, 2)(1.0, 1.0)
  S([1, 2])([1.0, 1.0])

  R([1, 2])(3.14)
  R(1, 2)(3.14)


  # all combinations
  S(1, 2, 3)(0.5, 0.5, 0.5)(CUBE(1))
  S([1, 2, 3])([0.5, 0.5, 0.5])(CUBE(1))

  T(1, 2, 3)(0.5, 0.5, 0.5)(CUBE(1))
  T([1, 2, 3])([0.5, 0.5, 0.5])(CUBE(1))

  STRUCT(S(1, 2, 3)(0.5, 0.5, 0.5), CUBE(1))
  STRUCT(S([1, 2, 3])([0.5, 0.5, 0.5]), CUBE(1))
  STRUCT(T(1, 2, 3)(0.5, 0.5, 0.5), CUBE(1))
  STRUCT(T([1, 2, 3])([0.5, 0.5, 0.5]), CUBE(1))

end


# ////////////////////////////////////////////////////////
function TestFenvs()

  begin

    @assert C(v -> sum(v))(1)(2) == 3
    @assert CAT([[1, 2], [3, 4]]) == [1, 2, 3, 4]
    @assert INV([[1.0, 0.0], [0.0, 1.0]]) == MatrixNd([[1.0, 0.0], [0.0, 1.0]])
    @assert EVERY(x -> x >= 0, [1, 2, 3, 4]) == true
    @assert EVERY(x -> x > 0, [1, -2, 3, 4]) == false
    @assert AND([true, true]) == true
    @assert AND([true, false]) == false
    @assert ID(10) == 10
    @assert DISTL([1.0, [2, 3, 4]]) == [(1.0, 2), (1.0, 3), (1.0, 4)]
    @assert DISTR([[1, 2, 3], 1.0]) == [(1, 1.0), (2, 1.0), (3, 1.0)]

    # preserve type
    @assert string(DISTL([100, [10.0, 20.0, 30.0]])) == "Any[(100, 10.0), (100, 20.0), (100, 30.0)]"
    @assert string(DISTR([[10.0, 20.0, 30.0], 100])) == "Any[(10.0, 100), (20.0, 100), (30.0, 100)]"
    @assert string(DISTL([100, 200, 300], [1.0, 2.0, 3.0])) == "Any[([100, 200, 300], 1.0), ([100, 200, 300], 2.0), ([100, 200, 300], 3.0)]"
    @assert string(DISTR([100, 200, 300], [0.1, 0.2, 0.3])) == "Any[(100, [0.1, 0.2, 0.3]), (200, [0.1, 0.2, 0.3]), (300, [0.1, 0.2, 0.3])]"

    @assert COMP([v -> [v; 3], v -> [v; 2], v -> [v; 1]])([0]) == [0, 1, 2, 3]
    @assert AA(x -> x * 2)([1, 2, 3]) == [2, 4, 6]
    @assert EQ([1, 1, 1]) == true
    @assert EQ([1, 1, 2]) == false
    @assert NEQ([1, 1, 2]) == true
    @assert NEQ([1, 1, 1]) == false
    @assert FILTER(LE(0))([-1, 0, 1, 2, 3, 4]) == [-1, 0]
    @assert FILTER(GE(0))([-1, 0, 1, 2, 3, 4]) == [0, 1, 2, 3, 4]
    @assert APPLY([x -> x * 2, 2]) == 4
    @assert INSL(x -> x[1] - x[2])([1, 2, 3]) == -4 # (1-2)-4
    @assert CONS([x -> x + 1, x -> x + 2])(0) == [1, 2]
    @assert IF([x -> x > 0, K(10), K(20)])(+1) == 10
    @assert IF([x -> x > 0, K(10), K(20)])(-1) == 20
    @assert LIFT(ADD)([cos, sin])(pi / 2.0) == 1.0
    @assert RAISE(ADD)([1, 2]) == 3
    @assert RAISE(ADD)([cos, sin])(pi / 2) == 1.0
    @assert ISNUM(0.0) == true
    @assert ISFUN(x -> x) && ISFUN(abs)
    @assert !ISFUN(3)
    @assert ISSEQOF(x -> isa(x, Int))([1, 2, 3]) == true
    @assert ISSEQOF(x -> isa(x, Int))([1, 2, 3.0]) == false
    @assert ISMAT([[1.0, 2.0], [3.0, 4.0]]) && ISMAT(MatrixNd(3))
    @assert ISMAT([1, 2, 3, 4]) == false
    @assert VECTSUM([[10, 11, 12], [0, 1, 2], [1, 1, 1]]) == [11, 13, 15]
    @assert VECTDIFF([[10, 11, 12], [0, 1, 2], [1, 1, 1]]) == [9, 9, 9]
    @assert MEANPOINT([[0, 0, 0], [1, 1, 1], [2, 2, 2]]) == [1.0, 1.0, 1.0]
    @assert SUM([1, 2, 3]) == 6
    @assert SUM([[1, 2, 3], [4, 5, 6]]) == [5, 7, 9]
    @assert SUM([[[1, 2], [3, 4]], [[10, 20], [30, 40]], [[100, 200], [300, 400]]]) == [[111, 222], [333, 444]]
    @assert DIFF(2) == -2
    @assert DIFF([1, 2, 3]) == -4
    @assert DIFF([[1, 2, 3], [1, 2, 3]]) == [0, 0, 0]
    @assert PROD([1, 2, 3, 4]) == 24
    @assert PROD([[1, 2, 3], [4, 5, 6]]) == 32
    @assert DIV([10.0, 2.0, 5.0]) == 1
    @assert REVERSE([1, 2, 3]) == [3, 2, 1]
    @assert REVERSE([1]) == [1]
    @assert TRANS([[1.0, 2.0], [3.0, 4.0]]) == [[1.0, 3.0], [2.0, 4.0]]
    @assert AR([[1, 2, 3], 0]) == [1, 2, 3, 0]
    @assert AL([0, [1, 2, 3]]) == [0, 1, 2, 3]
    @assert PROGRESSIVESUM([1, 2, 3, 4]) == [1, 3, 6, 10]
    @assert INTSTO(5) == [1, 2, 3, 4, 5]
    @assert FROMTO([1, 4]) == [1, 2, 3, 4]
    @assert S1([1, 2, 3]) == 1
    @assert S2([1, 2, 3]) == 2
    @assert N(3)(10) == [10, 10, 10]
    @assert DIESIS(3)(10) == [10, 10, 10]
    @assert NN(3)([10]) == [10, 10, 10]
    @assert DOUBLE_DIESIS(3)([10]) == [10, 10, 10]
    @assert AS(SEL)([1, 2, 3])([10, 11, 12]) == [10, 11, 12]
    @assert AC(SEL)([1, 2, 3])([10, 11, [12, [13]]]) == 13
    @assert CHARSEQ("hello") == ['h', 'e', 'l', 'l', 'o']
    @assert STRING(CHARSEQ("hello")) == "hello"
    @assert RANGE([1, 3]) == [1, 2, 3]
    @assert RANGE([3, 1]) == [3, 2, 1]
    @assert SIGN(10) == 1
    @assert SIGN(-10) == -1
    @assert TREE(x -> x[1] >= x[end] ? x[1] : x[end])([1, 2, 3, 4, 3, 2, 1]) == 4
    @assert TREE(x -> x[1] >= x[end] ? x[1] : x[end])([1, 2, 3, 4, 3, 2]) == 4
    @assert MERGE((x, y) -> x > y)([[1, 3, 4, 5], [2, 4, 8]]) == [1, 2, 3, 4, 4, 5, 8]
    @assert CASE([[LT(0), K(-1)], [C(EQ)(0), K(0)], [GT(0), K(+1)]])(-10) == -1
    @assert CASE([[LT(0), K(-1)], [C(EQ)(0), K(0)], [GT(0), K(+1)]])(0) == 0
    @assert CASE([[LT(0), K(-1)], [C(EQ)(0), K(0)], [GT(0), K(+1)]])(10) == 1
    @assert length(PERMUTATIONS([1, 2, 3])) == 6
    @assert toList(IDNT(0)) == []
    @assert toList(IDNT(2)) == [[1.0, 0.0], [0.0, 1.0]]
    @assert abs(SPLIT_2PI(4)[3] - PI) < 1e-4
    @assert VECTPROD([[1, 0, 0], [0, 1, 0]]) == [0.0, 0.0, 1.0]
    @assert VECTPROD([[0, 1, 0], [0, 0, 1]]) == [1.0, 0.0, 0.0]
    @assert VECTPROD([[0, 0, 1], [1, 0, 0]]) == [0.0, 1.0, 0.0]
    @assert VECTNORM([1, 0, 0]) == 1.0
    @assert INNERPROD([[1, 2, 3], [4, 5, 6]]) == 32.0
    @assert SCALARVECTPROD([2, [0, 1, 2]]) == [0, 2, 4]
    @assert SCALARVECTPROD([[0, 1, 2], 2]) == [0, 2, 4]
    @assert MIXEDPROD([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]) == 1.0
    @assert UNITVECT([2, 0, 0]) == [1.0, 0.0, 0.0]
    @assert UNITVECT([1, 1, 1]) == UNITVECT([2, 2, 2])
    @assert DIRPROJECT([1.0, 0.0, 0.0])([2.0, 0.0, 0.0]) == [2.0, 0.0, 0.0]
    @assert DIRPROJECT([1.0, 0.0, 0.0])([0.0, 1.0, 0.0]) == [0.0, 0.0, 0.0]
    @assert ORTHOPROJECT([1.0, 0.0, 0.0])([1.0, 1.0, 0.0]) == [0.0, 1.0, 0.0]
    @assert FACT(4) == 24
    @assert FACT(0) == 1
    @assert CHOOSE([7, 3]) == 35.0
    @assert TRACE(MatrixNd([[5.0, 0.0], [0.0, 10.0]])) == 15.0
    @assert toList(MATHOM(MatrixNd([[1, 2], [3, 4]]))) == [[1.0, 0.0, 0.0], [0.0, 1.0, 2.0], [0.0, 3.0, 4.0]]
    @assert SCALARMATPROD([10.0, [[1, 2], [3, 4]]]) == [[10.0, 20.0], [30.0, 40.0]]
    @assert ORTHO([[1, 0], [0, 1]]) == [[1.0, 0.0], [0.0, 1.0]]
    @assert SKEW([[1, 0], [0, 1]]) == [[0.0, 0.0], [0.0, 0.0]]
    @assert length(CART([[1, 2, 3], ['a', 'b'], [10, 11]])) == 12
    @assert length(POWERSET([1, 2, 3])) == 8
    @assert ISPOL(Cube(2)) == true
    @assert box(GRID([1, -1, 1])) == BoxNd([0.0], [3.0])
    @assert box(GRID([-1, 1, -1, 1])) == BoxNd([1.0], [4.0])
    @assert box(INTERVALS(10.0)(8)) == BoxNd([0.0], [10.0])
    @assert box(CUBOID([1, 2, 3])) == BoxNd([0.0, 0.0, 0.0], [1.0, 2.0, 3.0])
    @assert box(HPCSIMPLEX(3)) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    @assert RN(Cube(2)) == 2
    @assert DIM(Cube(2)) == 2
    @assert box(MKPOL([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], [[1, 2, 3, 4]])) == BoxNd([0.0, 0.0], [1.0, 1.0])
    @assert box((TRANSLATE(3)(2)(Cube(2)))) == BoxNd([0.0, 0.0, 2.0], [1.0, 1.0, 2.0])
    @assert box((TRANSLATE([1, 3])([1, 2])(Cube(2)))) == BoxNd([1.0, 0.0, 2.0], [2.0, 1.0, 2.0])
    @assert box((SCALE(3)(2)(Cube(3)))) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 2.0])
    @assert box((SCALE([3, 1])([4, 2])(Cube(3)))) == BoxNd([0.0, 0.0, 0.0], [2.0, 1.0, 4.0])
    @assert fuzzyEqual(box((ROTATE([1, 2])(PI / 2)(Cube(2)))), BoxNd([-1.0, 0.0], [0.0, +1.0]))
    @assert box((MAT([[1, 0, 0], [1, 1, 0], [2, 0, 1]])(Cube(2)))) == BoxNd([1.0, 2.0], [2.0, 3.0])

    # if you need to test shearing, enable this pirce
    if false

      # in 2dim
      VIEW(SHEARING(1)(0.5)(SQUARE(1)))
      VIEW(SHEARING(2)(0.5)(SQUARE(1)))

      # in 3dim
      VIEW(SHEARING(1)(0.2, 0.3)(CUBE(1)))
      VIEW(SHEARING(2)(0.2, 0.3)(CUBE(1)))
      VIEW(SHEARING(3)(0.2, 0.3)(CUBE(1)))

      #  in 4 dim
      SHEARING(1)(0.2, 0.3, 0.4)(CUBOID([1, 1, 1, 1]))
      SHEARING(2)(0.2, 0.3, 0.4)(CUBOID([1, 1, 1, 1]))
      SHEARING(3)(0.2, 0.3, 0.4)(CUBOID([1, 1, 1, 1]))
      SHEARING(4)(0.2, 0.3, 0.4)(CUBOID([1, 1, 1, 1]))

    end

    @assert box((STRUCT([
      Cube(2),
      T([1, 2])([1.0, 1.0]),
      T([1, 2])([1.0, 1.0]),
      Cube(2),
      Cube(2, 1.0, 2.0)]))) == BoxNd([0.0, 0.0], [4.0, 4.0])

    @assert box((STRUCT([
      T([1, 2])([1, 1]),
      T([1, 2])([1, 1]),
      Cube(2),
      T([1, 2])([1, 1]),
      T([1, 2])([1, 1]),
      Cube(2),
      Cube(2, 1.0, 2.0)]))) == BoxNd([2.0, 2.0], [6.0, 6.0])

    @assert box(JOIN([Cube(2, 0.0, 1.0)])) == BoxNd([0.0, 0.0], [1.0, 1.0])
    @assert POWER([2, 2]) == 4.0
    @assert box((POWER([Cube(2), Cube(1)]))) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    @assert SIZE(1)(Cube(2)) == 1.0
    @assert SIZE([1, 3])(SCALE([1, 2, 3])([1, 2, 3])(Cube(3))) == [1.0, 3.0]
    @assert MIN(1)(Cube(2)) == 0.0
    @assert MIN([1, 3])(TRANSLATE([1, 2, 3])([10.0, 20.0, 30.0])(Cube(3))) == [10.0, 30.0]
    @assert MAX(1)(Cube(2)) == 1.0
    @assert MAX([1, 3])(TRANSLATE([1, 2, 3])([10.0, 20.0, 30.0])(Cube(3))) == [11.0, 31.0]
    @assert MED(1)(Cube(2)) == 0.5
    @assert MED([1, 3])(Cube(3)) == [0.5, 0.5]
    @assert box((ALIGN([3, MAX, MIN])([Cube(3), Cube(3)]))) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 2.0])
    @assert box((TOP([Cube(3), Cube(3)]))) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 2.0])
    @assert box((BOTTOM([Cube(3), Cube(3)]))) == BoxNd([0.0, 0.0, -1.0], [1.0, 1.0, 1.0])
    @assert box((LEFT([Cube(3), Cube(3)]))) == BoxNd([-1.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    @assert box((RIGHT([Cube(3), Cube(3)]))) == BoxNd([0.0, 0.0, 0.0], [2.0, 1.0, 1.0])
    @assert box((UP([Cube(3, 0.0, 1.0), Cube(3, 5.0, 6.0)]))) == BoxNd([0.0, 0.0, 0.0], [6.0, 2.0, 1.0])
    @assert box((DOWN([Cube(3, 0.0, 1.0), Cube(3, 5.0, 6.0)]))) == BoxNd([0.0, -1.0, 0.0], [6.0, 1.0, 1.0])

    obj = Translate(Cube(3), [1.0, 2.0, 3.0])
    @assert box(BOX([1, 3])(obj)) == BoxNd([1.0, 3.0], [2.0, 4.0])

    obj = Translate(Cube(3), [1.0, 2.0, 3.0])
    @assert box(BOX(3)(obj)) == BoxNd([3.0], [4.0])

    @assert box((MAP([S1, S2])(Cube(2)))) == BoxNd([0.0, 0.0], [1.0, 1.0])
    @assert box((MAP(ID)(Cube(2)))) == BoxNd([0.0, 0.0], [1.0, 1.0])
    @assert box(CIRCUMFERENCE(1)(8)) == BoxNd([-1.0, -1.0], [1.0, 1.0])
    @assert box((RING([0.5, 1])([8, 8]))) == BoxNd([-1.0, -1.0], [1.0, 1.0])
    @assert box(CIRCLE(1.0)([8, 8])) == BoxNd([-1.0, -1.0], [1.0, 1.0])
    @assert fuzzyEqual(box(CYLINDER([1.0, 2.0])(8)), BoxNd([-1.0, -1.0, 0.0], [+1.0, +1.0, 2.0]))
    @assert box(SPHERE(1.0)([8, 8])) == BoxNd([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])
    @assert fuzzyEqual(box(TORUS([1.0, 2.0])([8, 8])), BoxNd([-2.0, -2.0, -0.5], [+2.0, +2.0, +0.5]))
    @assert fuzzyEqual(box((CONE([1.0, 3.0])(16))), BoxNd([-1.0, -1.0, 0.0], [+1.0, +1.0, 3.0]))

    println("All assert ok")

  end


  VIEW(ICOSPHERE())

  VIEW(icosphere(2))


  # some new tests
  begin


    # OK
    TRANS([[10, 20, 30], [1, 2, 3]])

    # OK
    hpc1 = Q(1.0)
    INSL(Power)([hpc1, hpc1, hpc1])

    # OK
    INSR(Power)([hpc1, hpc1, hpc1])

    # OK
    a = DISTR([[1, 2, 3], 40.0])
    @assert(a == [(1, 40.0), (2, 40.0), (3, 40.0)])

    # OK
    a = DISTL([40.0, [1, 2, 3]])
    @assert(a == [(40.0, 1), (40.0, 2), (40.0, 3)])

    # OK
    Q(1)

    # OK
    Q(1.0)

    # OK
    QUOTE([1, 2, -1, 3.0])

    # OK
    QUOTE([1, 2, -1, 3.0])

    # K is different from ID (see definitions from fenvs)
    #   ID returns the argument
    #   K(arg1)(arg2) returns always arg1
    @assert(ID(Q(1.0)) isa Plasm.Hpc)
    @assert(K(Q(1.0))("whatever") isa Plasm.Hpc)

    # OK
    hpc1 = Q(1.0)
    hpc2 = Power(hpc1, hpc1)
    @assert(dim(hpc2) == 2)

    # OK
    hpc3 = Power(hpc2, hpc1)
    @assert(dim(hpc3) == 3)
    # VIEW(hpc3) 

  end

  # ok
  begin

    obj = Struct([
      Translate(HpcSimplex(1), [0.0, 0.0, 0.0]),
      Translate(HpcSimplex(2), [1.0, 0.0, 0.0]),
      Translate(HpcSimplex(3), [2.0, 0.0, 0.0]),
      Translate(Cube(1), [0.0, 1.0, 0.0]),
      Translate(Cube(2), [1.0, 1.0, 0.0]),
      Translate(Cube(3), [2.0, 1.0, 0.0]),
    ])
    VIEW(obj, title="Example")

    TestArguments()
    TestCube()
    TestFenvsMkPol()
    TestSphere()
    TestMapSphere()
    TestTorus()
    TestBezier()
    TestCoonsPatch()
    TestRuledSurface()
    TestProfileProdSurface()
    TestRotationalSurface()
    TestConicalSurface()
    TestCubicHermite()
    TestSchegel3d()
    TestHermiteSurface()
    TestPermutahedron()
    TestBilinarSurface()
    TestBiquadraticSurface()
    TestBezierSurface()
    TestThinSolid()
    TestCubicSpline()
    # TestCylindricalSurface() 
    TestBezierManifold()
    TestOffset()
    TestEllipse()
    TestBezierStripe()
    TestDisplayNubSpline()
    TestDisplayNurbsSpline()

  end

  # example of `frame`
  begin
    VIEW(
      STRUCT(
        LINE([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], line_color=RED, line_width=DEFAULT_LINE_WIDTH),
        LINE([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], line_color=GREEN, line_width=DEFAULT_LINE_WIDTH),
        LINE([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], line_color=BLUE, line_width=DEFAULT_LINE_WIDTH),
      ),
      properties=Properties(
        "show_axis" => false,
        "title" => "Frame"
      )
    )
  end


  # example 2d
  begin
    VIEW(STRUCT([
      CUBOID([1, 1]),
      FRAME2([0.0, 0.0], [1.5, 1.5]),
    ]))
  end

  # example 2d
  begin
    VIEW(STRUCT([
      CUBOID([1, 1, 1]),
      FRAME3([0.0, 0.0, 0.0], [1.5, 1.5, 1.5]),
    ])
    )
  end

  begin

    # example WITH simplicial conversion
    set_config("map-convert-to-simplicial", true)
    obj = CIRCLE(1)([4, 1])
    #print(obj)
    @assert(all([length(hull) == 3 for hull in obj.childs[1].childs[1].hulls]))
    VIEW(obj, title="WITH simplicial conversion")
    set_config("map-convert-to-simplicial", false)

    # example WITHOUT simplicial conversion (NOw THE DEFAULT)
    obj = CIRCLE(1)([4, 1])
    #print(obj)
    @assert(all([length(hull) == 4 for hull in obj.childs[1].childs[1].hulls]))
    VIEW(obj, title="WITHOUT simplicial conversiom")
  end

  VIEW(
    CUBOID([1, 1, 1]),

    # view properties
    properties=Properties(
      "title" => "view with properties",

      # 3d pipeline position
      "pos" => Point3d(0.5, 0.5, 3.0),
      "dir" => Point3d(0.0, 0.0, -1.0),
      "vup" => Point3d(1.0, 0.0, 0.0),
      "znear" => 0.1,
      "zfar" => 10.0,

      # perspective fov
      "fov" => DEFAULT_FOV,

      #triangles, show/hide lines
      "show_lines" => false,

      #viewer background color
      "background_color" => DEFAULT_BACKGROUND_COLOR,

      # perspective or ortho projection
      "use_ortho" => DEFAULT_USE_ORTHO
    ))

  # BROKEN in julia
  # TestMinkowski()

  # broken in Python
  # tests.TestPolar()
  # tests.TestSolidify()
  # tests.TestDiff()

  println("TestFenvs ok")

end


# ///////////////////////////////////////////////////////
function TestHpc()
  TestViewer()
  TestComputeNormal()
  TestGoodTet()
  TestBox()
  TestHpcMat()
  TestHpcMkPol()
  TestHpcInternal()
  TestProperties()
  TestProperties2()
  Temple()
  Manhattan2D()
  Manhattan3D()
  TestFenvs()
  println("TestHpc ok")
end

TestHpc()



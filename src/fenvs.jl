using LinearAlgebra, Combinatorics

export SQRT, PI, COS, LEN, AND, OR, ToFloat64, C, ATAN2, MOD, ADD, MEANPOINT, SKEW,
  CAT, ISMAT, INV, EVERY, ID, K, DISTL, DISTR, COMP, AA, EQ, NEQ, LT, LE, GT, GE,
  ISGT, ISLT, ISGE, ISLE, BIGGER, SMALLER, FILTER, APPLY, INSR, INSL, BIGGEST, SMALLEST, CHARSEQ, STRING,
  CONS, IF, LIFT, RAISE, ISNUM, ISFUN, ISREAL, ISSEQ, ISSEQOF, VECTSUM, VECTDIFF, SUM,
  DIFF, PROD, SQR, DIV, REVERSE, TRANS,
  FIRST, LAST, TAIL, RTAIL, AR, AL, LIST, NOT, PROGRESSIVESUM,
  INTSTO, FROMTO, SEL, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, N, NN, DIESIS,
  DOUBLE_DIESIS, AS, AC, RANGE, SIGN, PRINT,
  PRINTPOL, TREE, MERGE, CASE, PERMUTATIONS, IDNT, SPLIT_2PI, VECTPROD,
  VECTNORM, INNERPROD, SCALARVECTPROD, MIXEDPROD, UNITVECT, DIRPROJECT, ORTHOPROJECT, FACT,
  ISREALVECT, ISFUNVECT, ISVECT, ISPOINT, CHOOSE, TRACE, MATHOM, SCALARMATPROD, MATDOTPROD, ORTHO, SUBSEQ,
  VECT2MAT, VECT2DTOANGLE, CART, POWERSET, ARC, ISPOL, PRINTPOL, PRINT, VIEW,
  GRID, QUOTE, INTERVALS, CUBOID, CUBE, HEXAHEDRON, SIMPLEX, RN, DIM, ISPOLDIM, MKPOL, MKPOLS,
  MK, UKPOL, UK, OPTIMIZE, TRANSLATE, T, SCALE, S, ROTATE, R, SHEARING, H,
  MAT, EMBED, STRUCT, HOLLOWCYL, SOLIDCYL,
  UNION, INTERSECTION, DIFFERENCE, XOR, CONVEXHULL,
  JOIN, POWER, SIZE, MIN, MAX, MED, ALIGN, TOP, BOTTOM, LEFT, RIGHT, UP, DOWN, BOX, MAP,
  CIRCLE_POINTS, CIRCUMFERENCE, NGON, RING, TUBE, CIRCLE, CYLINDER, CONE, TRUNCONE, DODECAHEDRON, ICOSAHEDRON, TETRAHEDRON,
  POLYPOINT, POLYLINE, TRIANGLESTRIPE, TRIANGLEFAN, MIRROR, POLYMARKER, BEZIER, BEZIERCURVE, COONSPATCH, RULEDSURFACE,
  PROFILEPRODSURFACE, ROTATIONALSURFACE, CYLINDRICALSURFACE, CONICALSURFACE, CUBICHERMITE, HERMITE, Q,
  EXTRUDE, MULTEXTRUDE, PROJECT, SPLITCELLS, EXTRACT_WIRES, SPLITPOLS, PERMUTAHEDRON, STAR, SCHLEGEL2D, SCHLEGEL3D,
  FINITECONE, PRISM, CROSSPOLYTOPE, OCTAHEDRON, ROTN, MKVERSORK, MKVECTOR, CUBICUBSPLINE, CUBICCARDINAL, SPLINE,
  JOINTS, BERNSTEINBASIS, TENSORPRODSURFACE, BILINEARSURFACE, BIQUADRATICSURFACE, HERMITESURFACE, BEZIERSURFACE,
  TENSORPRODSOLID, BEZIERMANIFOLD, LOCATE, RIF, FRACTALSIMPLEX, MESH, NU_GRID, SEGMENT, SOLIDIFY, EXTRUSION,
  EX, LEX, SEX, UKPOLF, POLAR, SWEEP, MINKOWSKI, OFFSET, THINSOLID, PLANE, RATIONALBEZIER, ELLIPSE, CURVE_NORMAL, DERBEZIER,
  BEZIERSTRIPE, BSPLINE, NUBSPLINE, DISPLAYNUBSPLINE, RATIONALBSPLINE, NURBSPLINE, DISPLAYNURBSPLINE, HOMO, PROPERTIES, SQUARE, LINE, MKPOINTS, FRAME2, FRAME3,
  COLOR, ICOSPHERE, icosphere, GRID1,
  TORUS, RING, SPHERE,
  IsPolytope, IsSimplex


import Base.sqrt
SQRT(f::Function) = x -> f(x)^(1 / 2)



import Base.-
-(f::Function, g::Function) = (x...) -> f(x...) - g(x...)

import Base.+
+(f::Function, g::Function) = (x...) -> f(x...) + g(x...)

import Base./
/(f::Function, g::Function) = (x...) -> f(x...) / g(x...)

import Base.*
*(f::Function, g::Function) = (x...) -> f(x...) * g(x...)

import Base.*
*(pol1::Hpc, pol2::Hpc) = Power(pol1, pol2)

import Base.^
^(f1::Function, f2::Function) = (x, y) -> f1(x)^f2(y)



PI = pi
COS = cos
LEN = length
AND = all
OR = any



# /////////////////////////////////////////////////////////////////
function ToFloat64(value)
  if isa(value, Vector)
    return [ToFloat64(it) for it in value]
  else
    return Float64(value)
  end
end

# /////////////////////////////////////////////////////////////////
function C(fun)
  function C1(arg1)
    function C2(arg2)
      return fun([arg1, arg2])
    end
    return C2
  end
  return C1
end

# /////////////////////////////////////////////////////////////////
ATAN2(l) = atan(l[2], l[1])

# /////////////////////////////////////////////////////////////////
MOD(l) = float(l[1] % l[2])

# /////////////////////////////////////////////////////////////////
function CAT(args)
  return reduce(vcat, args)
end

function ISMAT(v::Vector{Vector{Float64}})
  return true
end

function ISMAT(v::MatrixNd)
  return true
end

function ISMAT(v::Any)
  return false
end

# /////////////////////////////////////////////////////////////////
function INV(T::MatrixNd)
  return invert(T)
end

function INV(T::Vector{Vector{Float64}})
  return invert(MatrixNd(T))
end

# /////////////////////////////////////////////////////////////////
function EVERY(predicate, iterable)
  for x in iterable
    if !predicate(x)
      return false
    end
  end
  return true
end

# /////////////////////////////////////////////////////////////////
#function CURRY(fn, cargs..., ckwargs...)
#	function call_fn(fargs..., fkwargs...)
#		d = Dict(ckwargs...)
#		d = merge(d, fkwargs)
#		return fn(cargs..., fargs...; d...)
#	end
#	return call_fn
#end

# /////////////////////////////////////////////////////////////////
#function AND(list)
#	for i in list
#		if !i
#			return false
#		end
#	end
##	return true
#end

# /////////////////////////////////////////////////////////////////
ID(anyValue) = anyValue

# /////////////////////////////////////////////////////////////////
function K(AnyValue)
  function K0(obj)
    return AnyValue
  end
  return K0
end
TT = K(true)

# /////////////////////////////////////////////////////////////////
function DISTL(a, B::Vector)
  # see https://discourse.julialang.org/t/is-there-a-way-to-block-type-conversion/106302/6
  ret = []
  for b in B
    push!(ret, (a, b))
  end
  return ret
end

function DISTL(args)
  return DISTL(args...)
end

# /////////////////////////////////////////////////////////////////
function DISTR(A::Vector, b)
  # see https://discourse.julialang.org/t/is-there-a-way-to-block-type-conversion/106302/6
  ret = []
  for a in A
    push!(ret, (a, b))
  end
  return ret
end

function DISTR(args)
  return DISTR(args...)
end

# /////////////////////////////////////////////////////////////////
function COMP(Funs::Vector)
  function compose(f, g)
    function h(x)
      return f(g(x))
    end
    return h
  end
  return reduce(compose, Funs)
end


function COMP(a, b...)
  b = collect(b)
  v = Vector([a; b])
  return COMP(v)
end


# /////////////////////////////////////////////////////////////////
function AA(f)
  function AA0(args)
    return map(f, args)
  end
  return AA0
end

Eq(x, y) = x == y

# /////////////////////////////////////////////////////////////////
function EQ(List)
  for i in List
    if !(i == List[1])
      return false
    end
  end
  return true
end

# /////////////////////////////////////////////////////////////////
function NEQ(List)
  return !EQ(List)
end

LT(a) = function (b)
  return b < a
end

LE(a) = function (b)
  return b <= a
end

GT(a) = function (b)
  return b > a
end

GE(a) = function (b)
  return b >= a
end

# /////////////////////////////////////////////////////////////////
function ISGT(args)
  A, B = args
  return GT(A)(B)
end

function ISLT(args)
  A, B = args
  return LT(A)(B)
end

function ISGE(args)
  A, B = args
  return GE(A)(B)
end

function ISLE(args)
  A, B = args
  return LE(A)(B)
end

function BIGGER(args)
  A, B = args
  return A >= B ? A : B
end



function SMALLER(args)
  A, B = args
  return A <= B ? A : B
end




# /////////////////////////////////////////////////////////////////
function FILTER(predicate)
  function FILTER0(sequence)
    ret = []
    for item in sequence
      if predicate(item)
        push!(ret, item)
      end
    end
    return ret
  end
  return FILTER0
end

# /////////////////////////////////////////////////////////////////
function APPLY(args)
  f, x = args
  return f(x)
end

# /////////////////////////////////////////////////////////////////
function INSR(f)
  function INSR0(seq)
    N = length(seq)
    res = seq[N]
    for I in N-2:-1:1
      res = f([seq[I], res])
    end
    return res
  end
  return INSR0
end

SMALLEST = INSR(SMALLER)

# /////////////////////////////////////////////////////////////////
function INSL(f)
  function INSL0(seq)
    N = length(seq)
    res = seq[1]
    for I in 2:N
      res = f([res, seq[I]])
    end
    return res
  end
  return INSL0
end

BIGGEST = INSL(BIGGER)

# /////////////////////////////////////////////////////////////////
function CONS(Funs)
  return function (x)
    return [f(x) for f in Funs]
  end
end

# /////////////////////////////////////////////////////////////////
function IF(funs)
  function IF1(arg)
    f1, f2, f3 = funs
    return f1(arg) ? f2(arg) : f3(arg)
  end
  return IF1
end

# /////////////////////////////////////////////////////////////////
function LIFT(f)
  return function (funs)
    return COMP([f, CONS(funs)])
  end
end

# /////////////////////////////////////////////////////////////////
function RAISE(f)
  function RAISE0(args)
    return IF([ISSEQOF(ISFUN), LIFT(f), f])(args)
  end
  return RAISE0
end

# /////////////////////////////////////////////////////////////////
ISNUM(x) = isa(x, Int) || isa(x, Float64) || isa(x, Complex)

ISFUN(x) = isa(x, Function)
ISREAL(x) = isa(x, Float64)
ISSEQ(x) = isa(x, Array)
ISINT(x) = isa(x, Int)
ISINTPOS(x) = ISINT(x) && x > 0

# /////////////////////////////////////////////////////////////////
function ISSEQOF(type_checker)
  function ISSEQOF0(arg)
    if !isa(arg, Array)
      return false
    end
    for item in arg
      if !type_checker(item)
        return false
      end
    end
    return true
  end
  return ISSEQOF0
end

# /////////////////////////////////////////////////////////////////
function VECTSUM(vects)
  return [sum(x) for x in zip(vects...)]
end

function VECTDIFF(vects)
  return [l[1] - sum(l[2:end]) for l in zip(vects...)]
end

# /////////////////////////////////////////////////////////////////
function MEANPOINT(points)
  coeff = 1.0 / length(points)
  return [coeff * x for x in VECTSUM(points)]
end

# /////////////////////////////////////////////////////////////////
function SUM(args)
  if isa(args, Array) && ISPOL(args[1])
    return UNION(args)
  end
  if isa(args, Array) && ISNUM(args[1])
    return sum(args)
  end
  if isa(args, Array) && isa(args[1], Array)
    if isa(args[1][1], Array)
      return AA(VECTSUM)(zip(args...))
    else
      return VECTSUM(args)
    end
  end
  error("'+\' function has been applied to $args!")
end
ADD = SUM

# /////////////////////////////////////////////////////////////////
function DIFF(args)
  if isa(args, Array) && ISPOL(args[1])
    return DIFFERENCE(args)
  end
  if ISNUM(args)
    return -1 * args
  end
  if isa(args, Array) && ISNUM(args[1])
    return reduce((x, y) -> x - y, args)
  end
  if isa(args, Array) && isa(args[1], Array)
    if isa(args[1][1], Array)
      return AA(VECTDIFF)(zip(args...))
    else
      return VECTDIFF(args)
    end
  end
  error("\'-\' function has been applied to $args!")
end

# /////////////////////////////////////////////////////////////////
function PROD(args)
  if isa(args, Array) && ISPOL(args[1])
    return POWER(args)
  end
  if isa(args, Array) && ISSEQOF(ISNUM)(args)
    return reduce((x, y) -> x * y, args)
  end
  if isa(args, Array) && length(args) == 2 && ISSEQOF(ISNUM)(args[1]) && ISSEQOF(ISNUM)(args[2])
    return sum([a * b for (a, b) in zip(args[1], args[2])])
  end
  error("PROD function has been applied to $args!")
end

# /////////////////////////////////////////////////////////////////
SQR = RAISE(RAISE(PROD))([ID, ID])

function DIV(args)
  return reduce((x, y) -> x / Float64(y), args)
end

REVERSE(List) = reverse(List)

# /////////////////////////////////////////////////////////////////
function TRANS(List)
  if isa(List, MatrixNd)
    return List'
  elseif isa(List, Tuple) || isa(List, Array)
    ret = zip(List...)
    return [[jt for jt in it] for it in ret]
  else
    error("invalid argument")
  end
end

# /////////////////////////////////////////////////////////////////
FIRST(List) = List[1]
LAST(List) = List[end]
TAIL(List) = List[2:end]
RTAIL(List) = List[1:end-1]

# /////////////////////////////////////////////////////////////////
function AR(args)
  v, value = args
  return [v; value]
end

function AL(args)
  value, v = args
  return [value; v]
end

LIST(x) = [x]

function Not(x)
  return !x
end

NOT = AA(Not)

# /////////////////////////////////////////////////////////////////
function PROGRESSIVESUM(arg)
  ret, acc = [], 0
  for value in arg
    acc += value
    push!(ret, acc)
  end
  return ret
end

# /////////////////////////////////////////////////////////////////
function INTSTO(n)
  return collect(1:n)
end

function FROMTO(args)
  return collect(args[1]:args[2])
end

# /////////////////////////////////////////////////////////////////
SEL(n) = function (v)
  return v[n]
end

# /////////////////////////////////////////////////////////////////
S1 = SEL(1)
S2 = SEL(2)
S3 = SEL(3)
S4 = SEL(4)
S5 = SEL(5)
S6 = SEL(6)
S7 = SEL(7)
S8 = SEL(8)
S9 = SEL(9)
S10 = SEL(10)

# /////////////////////////////////////////////////////////////////
function N(n)
  return function (List)
    return [List for _ in 1:n]
  end
end

# /////////////////////////////////////////////////////////////////
function DIESIS(n)
  return function (List)
    return [List for _ in 1:n]
  end
end

# scrgiorgio: very dangerous to use D as a function name
# there is a risk of confusion between the function and the argument D
# D = DIESIS

# /////////////////////////////////////////////////////////////////
function NN(n)
  return function (v)
    return repeat(v, n)
  end
end

# /////////////////////////////////////////////////////////////////
function DOUBLE_DIESIS(n)
  return function (v)
    return repeat(v, n)
  end
end

DD = DOUBLE_DIESIS

# /////////////////////////////////////////////////////////////////
function AS(fun)
  return function (args)
    return COMP([CONS, AA(fun)])(args)
  end
end

# /////////////////////////////////////////////////////////////////
function AC(fun)
  return function (args)
    return COMP(AA(fun)(args))
  end
end

# /////////////////////////////////////////////////////////////////
function CHARSEQ(string)
  return collect(string)
end

STRING(charseq) = join(charseq)

# /////////////////////////////////////////////////////////////////
function RANGE(Pair)
  if Pair[end] - Pair[1] >= 0
    return collect(Pair[1]:Pair[end])
  end
  return collect(Pair[1]:-1:Pair[end])
end

# /////////////////////////////////////////////////////////////////
SIGN(Number) = Number >= 0 ? 1 : -1

# /////////////////////////////////////////////////////////////////
function PRINT(AnyValue)
  println(AnyValue)
  return AnyValue
end

function PRINTPOL(PolValue)
  println(PolValue)
  flush(stdout)
  return PolValue
end

# /////////////////////////////////////////////////////////////////
function TREE(f)
  function TREE1(v)
    N = length(v)
    if N == 1
      return v[1]
    end
    middle = trunc(Int, N / 2)
    return f([TREE1(v[1:middle]); TREE1(v[middle+1:N])])
  end
  return TREE1
end

# /////////////////////////////////////////////////////////////////
function MERGE(f)
  function MERGE1(v)
    a, b = v
    if length(a) == 0
      return b
    end
    if length(b) == 0
      return a
    end
    res = f(a[1], b[1])
    if !res
      return [a[1]; MERGE1([a[2:end], b])]
    else
      return [b[1]; MERGE1([a, b[2:end]])]
    end
  end
  return MERGE1
end

# /////////////////////////////////////////////////////////////////
function CASE(ListPredFuns)
  function CASE_NO_CURRIED(ListPredFuns, x)
    for p in ListPredFuns
      if p[1](x)
        return p[2](x)
      end
    end
  end
  return function (arg)
    return CASE_NO_CURRIED(ListPredFuns, arg)
  end
end

# /////////////////////////////////////////////////////////////////
function PERMUTATIONS(SEQ)
  if length(SEQ) <= 1
    return [SEQ]
  end
  ret = []
  for i in 1:length(SEQ)
    element = SEQ[i]
    rest = PERMUTATIONS([SEQ[1:i-1]; SEQ[i+1:end]])
    for r in rest
      push!(ret, [element; r])
    end
  end
  return ret
end

# /////////////////////////////////////////////////////////////////
function IDNT(N)
  return MatrixNd(N)
end

# /////////////////////////////////////////////////////////////////
function SPLIT_2PI(N)
  delta = 2 * PI / N
  return [i * delta for i in 0:N-1]
end

# /////////////////////////////////////////////////////////////////
function VECTPROD(args)
  a, b = args
  ax, ay, az = [Float64(it) for it in a]
  bx, by, bz = [Float64(it) for it in b]
  return [
    ay * bz - by * az,
    az * bx - bz * ax,
    ax * by - bx * ay
  ]
end

# /////////////////////////////////////////////////////////////////
function VECTNORM(u)
  return sqrt(sum(x * x for x in u))
end

# /////////////////////////////////////////////////////////////////
function INNERPROD(args)
  return COMP([COMP([RAISE(SUM), AA(RAISE(PROD))]), TRANS])(args)
end

# /////////////////////////////////////////////////////////////////
function SCALARVECTPROD(args)
  s, l = args
  if !isa(l, Array)
    s, l = l, s
  end
  return [s * l[i] for i in 1:length(l)]
end

# /////////////////////////////////////////////////////////////////
function MIXEDPROD(args)
  A, B, C = args
  return INNERPROD([VECTPROD([A, B]), C])
end

# /////////////////////////////////////////////////////////////////
function UNITVECT(V)
  norm = VECTNORM(V)
  return [coord / norm for coord in V]
end

# /////////////////////////////////////////////////////////////////
function DIRPROJECT(E)
  function DIRPROJECT0(V)
    return SCALARVECTPROD([(INNERPROD([E, V])), E])
  end
  return DIRPROJECT0
end

# /////////////////////////////////////////////////////////////////
function ORTHOPROJECT(E)
  function ORTHOPROJECT0(V)
    return VECTDIFF([V, DIRPROJECT((E))(V)])
  end
  return ORTHOPROJECT0
end


# /////////////////////////////////////////////////////////////////
FACT(n) = n > 0 ? *(1:big(n)...) : 1

ISREALVECT = ISSEQOF(ISREAL)
ISFUNVECT = ISSEQOF(ISFUN)
ISVECT = COMP([OR, CONS([ISREALVECT, ISFUNVECT])])
ISPOINT = ISVECT

# /////////////////////////////////////////////////////////////////
function CHOOSE(args)
  N, K = args
  return FACT(N) / float(FACT(K) * FACT(N - K))
end

# /////////////////////////////////////////////////////////////////
function TRACE(MATRIX)
  acc = 0.0
  N = dim(MATRIX)
  for I in 1:N
    acc += MATRIX[I, I]
  end
  return acc
end

# /////////////////////////////////////////////////////////////////
function MATHOM(T)
  N = dim(T)
  ret = MatrixNd(N + 1)
  ret[2:N+1, 2:N+1] = T[1:N, 1:N]
  return ret
end

# /////////////////////////////////////////////////////////////////
function SCALARMATPROD(args)
  scalar, mat = float(args[1]), args[2]
  return [[scalar * coord for coord in row] for row in mat]
end

# /////////////////////////////////////////////////////////////////
function MATDOTPROD(args)
  return COMP([INNERPROD, AA(CAT)])(args)
end

# /////////////////////////////////////////////////////////////////
function ORTHO(matrix)
  return SCALARMATPROD([0.5, SUM([matrix, TRANS(matrix)])])
end

# /////////////////////////////////////////////////////////////////
function SKEW(matrix)
  return SCALARMATPROD([0.5, DIFF([matrix, TRANS(matrix)])])
end

# /////////////////////////////////////////////////////////////////
function SUBSEQ(I_J)
  function SUBSEQ0(SEQ)
    return SEQ[I_J[1]:I_J[2]]
  end
  return SUBSEQ0
end

# /////////////////////////////////////////////////////////////////
function VECT2MAT(v)
  n = length(v)
  return [[r != c ? 0 : v[r] for c in 1:n] for r in 1:n]
end

# /////////////////////////////////////////////////////////////////
function VECT2DTOANGLE(v)
  v = UNITVECT(v)
  return acos(v[1]) * (v[2] >= 0 ? 1 : -1)
end

# /////////////////////////////////////////////////////////////////
function CART(l)
  return vec(collect(Iterators.product(l...)))
end

# /////////////////////////////////////////////////////////////////
function POWERSET(l)
  return collect(powerset([1, 2, 3]))
end

# /////////////////////////////////////////////////////////////////
function ARC(args)
  degrees, cents = args
  return PI * (degrees + cents) / (100.0 * 180.0)
end

# /////////////////////////////////////////////////////////////////
function ISPOL(obj)
  return isa(obj, Hpc)
end




# /////////////////////////////////////////////////////////////////
VIEW = View

function LINE(p1, p2; line_color=Point4d(1.0, 1.0, 1.0, 1.0), line_width=1)
  return PROPERTIES(
    MKPOL([p1, p2], [[1, 2]]),
    Properties(
      "line_color" => line_color,
      "line_width" => line_width
    )
  )
end

function MKPOINTS(points)
  return MKPOL(points, [[I] for I in 1:length(points)])
end

function FRAME2(p1=[0.0, 0.0], p2=[1.0, 1.0])
  return STRUCT(
    LINE([p1[1], p1[2]], [p2[1], p1[2]], line_color=RED, line_width=2),
    LINE([p1[1], p1[2]], [p1[1], p2[2]], line_color=GREEN, line_width=2)
  )
end

function FRAME3(p1=[0.0, 0.0, 0.0], p2=[1.0, 1.0, 1.0])
  return STRUCT(
    LINE([p1[1], p1[2], p1[3]], [p2[1], p1[2], p1[3]], line_color=RED, line_width=2),
    LINE([p1[1], p1[2], p1[3]], [p1[1], p2[2], p1[3]], line_color=GREEN, line_width=2),
    LINE([p1[1], p1[2], p1[3]], [p1[1], p1[2], p2[3]], line_color=BLUE, line_width=2),
  )
end






# /////////////////////////////////////////////////////////////////
function GRID(sequence)
  cursor = 0
  points = [[0.0]]
  hulls = Vector{Vector{Int}}()
  N = 1
  for value in sequence
    cursor += abs(value)
    push!(points, [cursor])
    N += 1
    if value >= 0
      push!(hulls, [N - 1, N])
    end
  end
  return MkPol(points, hulls)
end
QUOTE = GRID

# Q = COMP([QUOTE, IF([ISSEQ, ID, CONS([ID])])])

# /////////////////////////////////////////////////////////////////
function INTERVALS(A)
  A = Float64(A)
  function INTERVALS0(N::Int64)
    return QUOTE([A / N for I in 1:N])
  end
  return INTERVALS0
end

# /////////////////////////////////////////////////////////////////
function CUBOID(vs)
  return Scale(Cube(length(vs)), [Float64(it) for it in vs])
end

function CUBE(size)
  return CUBOID([Float64(size), Float64(size), Float64(size)])
end

SQUARE(d) = CUBOID([d, d])


# /////////////////////////////////////////////////////////////////
function HEXAHEDRON()
  return Cube(3, -1.0 / sqrt(3.0), +1.0 / sqrt(3.0))
end

# /////////////////////////////////////////////////////////////////
"""
    SIMPLEX(dim::Int)::Hpc
Generator of the standard `Hpc` simplex object of dimension `dim`.

Simplex object of `Hpc` type with arbitrary dimension `dim ≥ 1`. 
It is the convex combination of ``d+1`` affinely independent points.

# Examples
```
julia> SIMPLEX(1)
Hpc(MatrixNd(2), Geometry[Geometry([[0.0], [1.0]], hulls=[[1, 2]])])

julia> SIMPLEX(2)
Hpc(MatrixNd(3), Geometry[Geometry([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], hulls=[[1, 2, 3]])])

julia> SIMPLEX(3)
Hpc(MatrixNd(4), Geometry[Geometry([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], hulls=[[1, 2, 3, 4]])])
```
"""
function SIMPLEX(dim)
  return Simplex(dim)
end

RN(pol) = dim(pol)

DIM(pol) = dim(pol)

# /////////////////////////////////////////////////////////////////
function ISPOLDIM(dims)
  function ISPOLDIM1(pol)
    d = dims[1]
    n = dims[2]
    return (d == DIM(pol)) && (n == RN(pol))
  end
  return ISPOLDIM1
end

# /////////////////////////////////////////////////////////////////
function MKPOL(points::Vector{Vector{Float64}}, hulls::Vector{Vector{Int}}, __pols=Nothing)
  return MkPol(points, hulls)
end

function MKPOL(V::Matrix{Float64}, hulls::Vector{Vector{Int}}, __pols=Nothing)
  W = [V[:, k] for k = 1:size(V)[2]]
  return MkPol(W, hulls)
end

# deprecated
function MKPOL(args)
  return MKPOL(args[1], args[2])
end

# /////////////////////////////////////////////////////////////////
MK = COMP([MKPOL, CONS([LIST, K([[1]]), K([[1]])])])

# /////////////////////////////////////////////////////////////////
function CONVEXHULL(points)
  return MKPOL(points, [collect(1:length(points))], [[1]])
end

# /////////////////////////////////////////////////////////////////
function MKPOLS(V::Vector{Vector{Float64}}, hulls::Vector{Vector{Int}})::Hpc
  out = STRUCT(AA(MKPOL)(DISTL(V, AA(LIST)(hulls))))
  return out
end

function MKPOLS(V::Matrix{Float64}, hulls::Vector{Vector{Int}})
  W = [V[:, k] for k = 1:size(V, 2)]
  return MKPOLS(W, hulls)
end

function MKPOLS(V::Union{Vector{Vector{Float64}},Matrix{Float64}}, cells::Dict{Symbol,AbstractArray})
  v = []
  for (__symbol, hulls) in cells
    push!(v, MKPOLS(V, hulls))
  end
  return STRUCT(v)
end

# /////////////////////////////////////////////////////////////////
export SKELETON
function SKELETON(ord::Int)
   function SKELETON0(pol::Hpc)
      larpol = LAR(pol)
      if ord == 1
         return MKPOLS(larpol.V, larpol.C[:EV])
      elseif ord == 2
         return MKPOLS(larpol.V, larpol.C[:FV])
      elseif ord == 3
         return MKPOLS(larpol.V, larpol.C[:CV])
      else
         error("not yet implemented")
      end
   end
   return SKELETON0
end

# /////////////////////////////////////////////////////////////////
function UKPOL(pol)
  points, hulls = UkPol(pol)
  ret = Vector{Any}()
  push!(ret, points)
  push!(ret, hulls)
  push!(ret, [[1]])
  return ret
end

UK = COMP([COMP([S1, S1]), UKPOL])

# /////////////////////////////////////////////////////////////////
OPTIMIZE(pol) = pol


# //////////////////////////////////////////////////////////////////////////////
"""
    IsPolytope
Plasm predicate `Expr -> Bool` in pure FL style.

Polytopes are the generalization of three-dimensional polyhedra to any number of dimensions.
"""
IsPolytope = AND ∘ CONS([  # pure FL style
  ISPOL,
  EQ ∘ CONS([LEN ∘ S2 ∘ UKPOL, K(1)])
])

# //////////////////////////////////////////////////////////////////////////////
"""
    IsSimplex
Plasm predicate `Expr -> Bool` in pure FL style.

generalization of the notion of a triangle or tetrahedron to arbitrary dimensions.
"""
IsSimplex = AND ∘ CONS([  # pure FL style
  IsPolytope,
  EQ ∘ CONS([LEN ∘ S1 ∘ UKPOL, RN + K(1)])
])


# /////////////////////////////////////////////////////////////////
function TRANSLATE(axis, values, pol::Hpc)
  axis = ISNUM(axis) ? [axis] : axis
  values = ISNUM(values) ? [values] : values
  vt = [0.0 for I in 1:max(axis...)]
  for (a, t) in zip(axis, values)
    vt[a] = t
  end
  return Translate(pol, vt)
end

function TRANSLATE(axis)
  function TRANSLATE1(values)
    function TRANSLATE2(pol::Hpc)
      return TRANSLATE(axis, values, pol)
    end
    return TRANSLATE2
  end
  return TRANSLATE1
end

function TRANSLATE(a, b...)
  b = collect(b)
  axis = [a; b]
  function TRANSLATE1(c, d...)
    d = collect(d)
    values = [c; d]
    function TRANSLATE2(pol::Hpc)
      return TRANSLATE(axis, values, pol)
    end
    return TRANSLATE2
  end
  return TRANSLATE1
end

T = TRANSLATE

# /////////////////////////////////////////////////////////////////
function SCALE(axis, values, pol::Hpc)
  axis = ISNUM(axis) ? [axis] : axis
  values = ISNUM(values) ? [values] : values
  vs = [1.0 for I in 1:max(axis...)]
  for (a, t) in zip(axis, values)
    vs[a] = t
  end
  return Scale(pol, vs)
end

function SCALE(axis)
  function SCALE1(values)
    function SCALE2(pol::Hpc)
      return SCALE(axis, values, pol)
    end
    return SCALE2
  end
  return SCALE1
end

function SCALE(a, b...)
  b = collect(b)
  axis = [a; b]
  function SCALE1(c, d...)
    d = collect(d)
    values = [c; d]
    function SCALE2(pol::Hpc)
      return SCALE(axis, values, pol)
    end
    return SCALE2
  end
  return SCALE1
end

S = SCALE

# /////////////////////////////////////////////////////////////////
function ROTATE(plane_indexes, angle, pol::Hpc)
  return Rotate(pol, plane_indexes[1], plane_indexes[2], angle)
end

function ROTATE(plane_indexes)
  function ROTATE1(angle)
    function ROTATE2(pol::Hpc)
      return ROTATE(plane_indexes, angle, pol)
    end
    return ROTATE2
  end
  return ROTATE1
end

function ROTATE(a, b...)
  b = collect(b)
  plane_indexes = [a; b]
  function ROTATE1(angle)
    function ROTATE2(pol::Hpc)
      return ROTATE(plane_indexes, angle, pol)
    end
    return ROTATE2
  end
  return ROTATE1
end


R = ROTATE

# /////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////
function SHEARING(column)
  function SHEARING1(a, b...)
    b = collect(b)
    shear_list = [a; b]
    function SHEARING2(pol::Hpc)
      return SHEAR(column, shear_list, pol)
    end
    return SHEARING2
  end
  return SHEARING1
end
H = SHEARING

function SHEAR(column, shear_list, pol::Hpc)
  values = ISNUM(shear_list) ? [shear_list] : shear_list
  vh = [0.0 for I in 1:length(values)+2]
  for a in 2:length(vh)-1
    vh[column+1] = 1
    if a < column + 1
      vh[a] = values[a-1]
    elseif a >= column + 1
      vh[a+1] = values[a-1]
    end
  end
  return Shear(pol, column, vh)
end

function Shear(column, vh)
  T = MatrixNd(length(vh))
  for I in 1:dim(T)
    T[I, column+1] = vh[I]
  end
  return T
end


function Shear(self::Hpc, column, vh::Vector{Float64})
  return Hpc(Shear(column, vh), [self])
end

H = SHEARING

# /////////////////////////////////////////////////////////////////
function MAT(matrix)
  function MAT0(pol)
    return Transform(pol, MatrixNd(matrix))
  end
  return MAT0
end

# /////////////////////////////////////////////////////////////////
# add homo coordinate
function HOMO(T)
  N = size(T, 1)
  ret = MatrixNd(N + 1)
  for I in 1:N
    for J in 1:N
      ret[I+1, J+1] = T[I, J]
    end
  end
  return ret
end

# /////////////////////////////////////////////////////////////////
function EMBED(up_dim)
  function EMBED1(pol)
    return Hpc(MatrixNd(dim(pol) + up_dim + 1), [pol])
  end
  return EMBED1
end

# /////////////////////////////////////////////////////////////////
function STRUCT(seq::Vector, nrec::Int=0)

  if isempty(seq)
    error("STRUCT must be applied to a not empty list")
  end

  if nrec == 0
    seq = copy(seq)
  end

  # collect all geometry wihout transformations
  pols = Vector{Hpc}()
  while !isempty(seq) && ISPOL(seq[1])
    push!(pols, seq[1])
    seq = seq[2:end]
  end

  # callect all transformation on the right
  transformations = []
  while !isempty(seq) && ISFUN(seq[1])
    push!(transformations, seq[1])
    seq = seq[2:end]
  end

  if !isempty(seq) && !ISPOL(seq[1]) && !ISFUN(seq[1])
    error("STRUCT arguments not valid, not all elements are pols or transformations")
  end

  # recursive on the right, apply transformations
  if !isempty(seq)
    @assert(ISPOL(seq[1]))
    child = STRUCT(seq, nrec + 1)
    @assert(ISPOL(child))
    if !isempty(transformations)
      child = COMP(transformations)(child)
    end
    push!(pols, child)
  end

  if isempty(pols)
    error("Cannot find geometry in STRUCT, found only transformations")
  end
  return Struct(pols)
end

# struct with a sequence
function STRUCT(a, b...)
  b = collect(b)
  v = Vector([a; b])
  return STRUCT(v)
end

# /////////////////////////////////////////////////////////////////////////////
struct BspFace
  vertices::Vector{Vector{Float64}}
  plane::Vector{Float64}

  function BspFace(vertices, plane=nothing)

    # automatic computation of plane
    if plane == nothing
      dim = length(vertices[1])
      if dim == 2
        @assert length(vertices) == 2
        x1, y1 = vertices[1][1], vertices[1][2]
        x2, y2 = vertices[2][1], vertices[2][2]
        normal = [(y2 - y1), -(x2 - x1)]
      elseif dim == 3
        normal = ComputeTriangleNormal(vertices[1], vertices[2], vertices[3])
      else
        error("todo")
      end

      w = sum([normal[i] * vertices[1][i] for i in 1:dim])
      plane = [normal; -w]
    end
    new(vertices, plane)
  end

end


function splitEdge(plane::Vector{Float64}, vi::Vector{Float64}, vj::Vector{Float64})
  dim = length(vi)
  valuev1 = abs(plane[end] + sum([plane[i] * vi[i] for i in 1:dim]))
  valuev2 = abs(plane[end] + sum([plane[i] * vj[i] for i in 1:dim]))
  alpha = 1.0 / (valuev1 + valuev2)
  beta = valuev1 * alpha
  alpha = alpha * valuev2
  return [alpha * vi[i] + beta * vj[i] for i in 1:dim]
end

function splitFace(self::BspFace, plane::Vector{Float64}, EPSILON::Float64=1e-5)
  dim = length(plane) - 1
  COPLANAR, ABOVE, BELOW, SPANNING = 0, 1, 2, 3
  ptype, types = COPLANAR, []
  for v in self.vertices
    @assert length(v) == dim
    t = plane[end] + sum([plane[i] * v[i] for i in 1:dim])
    if t < -EPSILON
      type = BELOW
    elseif t > EPSILON
      type = ABOVE
    else
      type = COPLANAR
    end
    ptype |= type
    push!(types, type)
  end
  if ptype == BELOW
    return [self, nothing, nothing, nothing]
  end
  if ptype == ABOVE
    return [nothing, nothing, nothing, self]
  end
  if ptype == COPLANAR
    if sum([plane[i] * self.plane[i] for i in 1:dim+1]) > 0
      return [nothing, nothing, self, nothing]
    else
      return [nothing, self, nothing, nothing]
    end
  end
  @assert ptype == SPANNING
  b, f = [], []
  if dim == 2
    @assert length(self.vertices) == 2
    ti, tj = types[1], types[2]
    vi, vj = self.vertices[1], self.vertices[2]
    if ti != BELOW
      push!(f, vi)
    end
    if ti != ABOVE
      push!(b, vi)
    end
    if tj != BELOW
      push!(f, vj)
    end
    if tj != ABOVE
      push!(b, vj)
    end
    if (ti | tj) == SPANNING
      v = splitEdge(plane, vi, vj)
      push!(b, v)
      push!(f, v)
    end
  elseif dim == 3
    for i in 1:length(self.vertices)
      j = (i + 1) % length(self.vertices)
      ti, tj = types[i], types[j]
      vi, vj = self.vertices[i], self.vertices[j]
      if ti != BELOW
        push!(f, vi)
      end
      if ti != ABOVE
        push!(b, vi)
      end
      if (ti | tj) == SPANNING
        v = splitEdge(plane, vi, vj)
        push!(b, v)
        push!(f, v)
      end
    end
  else
    error("not supported")
  end
  @assert length(b) >= dim && length(f) >= dim
  return [BspFace(b, self.plane), nothing, nothing, BspFace(f, self.plane)]
end

# //////////////////////////////////////////////////////////////////////////////////////
struct Bsp

  plane::Vector{Float64}
  faces::Vector{BspFace}
  below::Bsp
  above::Bsp

  function Bsp()
    new()
  end

end

function getFaces(self::Bsp)
  if self == nothing
    return []
  end
  return [self.faces; getFaces(self.below); getFaces(self.above)]
end

function insertFaces(self::Bsp, faces::Vector{BspFace})
  if length(faces) == 0
    return self
  end

  if self.plane == nothing
    @assert self.below === nothing && self.above === nothing
    self.plane = faces[1].plane
  end

  below, above = [], []
  for face in faces
    b, cb, ca, a = splitFace(face, self.plane)
    if b != nothing
      push!(below, b)
    end
    if cb != nothing
      push!(self.faces, cb)
    end
    if ca != nothing
      push!(self.faces, ca)
    end
    if a != nothing
      push!(above, a)
    end
  end

  if !isempty(above)
    if self.above === nothing
      self.above = Bsp()
    end
    insertFaces(self.above, above)
  end

  if !isempty(below)
    if self.below === nothing
      self.below = Bsp()
    end
    insertFaces(self.below, below)
  end
end

function clipFaces(self::Bsp, faces)
  if self.plane === nothing
    return faces
  end
  below, above = [], []
  for p in faces
    b, cb, ca, a = splitFace(p, self.plane)
    if b !== nothing
      push!(below, b)
    end
    if cb !== nothing
      push!(below, cb)
    end
    if ca !== nothing
      push!(above, ca)
    end
    if a !== nothing
      push!(above, a)
    end
  end
  below = self.below !== nothing ? clipFaces(self.below, below) : []
  above = self.above !== nothing ? clipFaces(self.above, above) : []
  return above + below
end

function clipTo(self::Bsp, other)
  self.faces = clipFaces(other, self.faces)
  if self.below !== nothing
    clipTo(self.below, other)
  end
  if self.above !== nothing
    clipTo(self.above, other)
  end
end

function Complement(bsp::Bsp)
  if bsp === nothing
    return nothing
  end
  ret = Bsp()
  if bsp.plane !== nothing
    ret.plane = [-1 * it for it in bsp.plane]
  end
  for p in bsp.faces
    new_p = BspFace(reverse(p.vertices), [-1 * c for c in p.plane])
    push!(ret.faces, new_p)
  end
  ret.below = Complement(bsp.above)
  ret.above = Complement(bsp.below)
  return ret
end

function Union(a::Bsp, b::Bsp)
  clipTo(a, b)
  clipTo(b, a)
  b = Bsp.Complement(b)
  clipTo(b, a)
  b = Bsp.Complement(b)
  insertFaces(a, getFaces(b))
  return a
end

function Intersection(a::Bsp, b::Bsp)
  return Bsp.Complement(Bsp.Union(Bsp.Complement(a), Bsp.Complement(b)))
end

function Difference(a::Bsp, b::Bsp)
  return Bsp.Complement(Bsp.Union(Bsp.Complement(a), b))
end

function Xor(a::Bsp, b::Bsp)
  return Bsp.Union(Bsp.Intersection(a, Bsp.Complement(b)), Bsp.Intersection(Bsp.Complement(a), b))
end

function fromHpc(hpc::Hpc)
  ret = Bsp()
  faces = []
  for (T, properties, obj) in toList(ToBoundaryForm(hpc))
    for hull in obj.hulls
      points = [transformPoint(T, obj.points[I]) for I in hull]
      push!(faces, BspFace(points))
    end
  end
  insertFaces(ret, faces)
  return ret
end

function toHpc(self::Bsp)
  batches, faces = [], getFaces(self)
  dim = self.plane !== nothing ? length(self.plane) - 1 : 0
  if dim == 0
    return Hpc()
  end
  @assert dim == 1 || dim == 2 || dim == 3
  points, hulls = [], []
  for face in faces
    if dim == 1
      @assert length(face.vertices) == 1
      N = length(points)
      points += face.vertices
      push!(hulls, collect(N:N+length(face.vertices)-1))
    elseif dim == 2
      @assert length(face.vertices) == 2
      N = length(points)
      points += face.vertices
      push!(hulls, collect(N:N+length(face.vertices)-1))
    elseif dim == 3
      @assert length(face.vertices) >= 3
      for I in 2:length(face.vertices)-1
        N = length(points)
        points += [face.vertices[1], face.vertices[I], face.vertices[I+1]]
        push!(hulls, collect(N:N+2))
      end
    else
      error("not supported")
    end
  end
  return MkPol(points, hulls)
end

function UNION(objs::Vector{Hpc})
  objs = [fromHpc(obj) for obj in objs]
  res = objs[1]
  for I in 2:length(objs)
    res = Union(res, objs[I])
  end
  return toHpc(res)
end

function INTERSECTION(objs::Vector{Hpc})
  objs = [fromHpc(obj) for obj in objs]
  res = objs[1]
  for I in 2:length(objs)
    res = Intersection(res, objs[I])
  end
  return toHpc(res)
end

function DIFFERENCE(objs::Vector{Hpc})
  objs = [fromHpc(obj) for obj in objs]
  res = objs[1]
  for I in 2:length(objs)
    res = Difference(res, objs[I])
  end
  return toHpc(res)
end

function XOR(objs::Vector{Hpc})
  objs = [fromHpc(obj) for obj in objs]
  res = objs[1]
  for I in 2:length(objs)
    res = Xor(res, objs[I])
  end
  return toHpc(res)
end

# ///////////////////////////////////////////////////////////
function JOIN(pol_list)
  if ISPOL(pol_list)
    pol_list = [pol_list]
  end
  return Join(pol_list)
end

# ///////////////////////////////////////////////////////////
function POWER(objs_list)
  if !isa(objs_list, Vector) || length(objs_list) != 2
    error("POWER can only be applied to a list of 2 arguments")
  end
  if ISNUM(objs_list[1]) && ISNUM(objs_list[2])
    return objs_list[1]^objs_list[2]
  end
  return Power(objs_list[1], objs_list[2])
end

# ///////////////////////////////////////////////////////////
function SIZE(sel)
  function SIZE1(pol)
    S = size(box(pol))
    return isa(sel, Vector) ? [S[i] for i in sel] : S[sel]
  end
  return SIZE1
end

# ///////////////////////////////////////////////////////////
function MIN(sel)
  function MIN1(pol)
    b = box(pol)
    return isa(sel, Vector) ? [b.p1[i] for i in sel] : b.p1[sel]
  end
  return MIN1
end

# ///////////////////////////////////////////////////////////
function MAX(sel)
  function MAX1(pol)
    b = box(pol)
    return isa(sel, Vector) ? [b.p2[i] for i in sel] : b.p2[sel]
  end
  return MAX1
end

# ///////////////////////////////////////////////////////////
function MED(sel)
  function MED1(pol)
    c = center(box(pol))
    return isa(sel, Vector) ? [c[i] for i in sel] : c[sel]
  end
  return MED1
end

# ///////////////////////////////////////////////////////////
function ALIGN(args)
  function ALIGN0(args, pols)
    pol1, pol2 = pols
    box1, box2 = box(pol1), box(pol2)
    if isa(args, Vector) && !isempty(args) && ISNUM(args[1])
      args = [args]
    end
    max_index = max([index for (index, pos1, po2) in args]...)
    vt = zeros(max_index)
    for (index, pos1, pos2) in args
      p1 = ifelse(pos1 == MIN, box1.p1, ifelse(pos1 == MAX, box1.p2, center(box1)))
      p1 = ifelse(index <= length(p1), p1[index], 0.0)
      p2 = ifelse(pos2 == MIN, box2.p1, ifelse(pos2 == MAX, box2.p2, center(box2)))
      p2 = ifelse(index <= length(p2), p2[index], 0.0)
      vt[index] -= (p2 - p1)
    end
    return Struct([pol1, Translate(pol2, vt)])
  end
  return pols -> ALIGN0(args, pols)
end
TOP = ALIGN([[3, MAX, MIN], [1, MED, MED], [2, MED, MED]])
BOTTOM = ALIGN([[3, MIN, MAX], [1, MED, MED], [2, MED, MED]])
LEFT = ALIGN([[1, MIN, MAX], [3, MIN, MIN]])
RIGHT = ALIGN([[1, MAX, MIN], [3, MIN, MIN]])
UP = ALIGN([[2, MAX, MIN], [3, MIN, MIN]])
DOWN = ALIGN([[2, MIN, MAX], [3, MIN, MIN]])

# ///////////////////////////////////////////////////////////
function BOX(sel)
  function BOX0(pol)
    if !isa(sel, Vector)
      sel = [sel]
    end
    dim = length(sel)
    b = box(pol)
    vt = [b.p1[i] for i in sel]
    vs = [size(b)[i] for i in sel]
    return Translate(Scale(Cube(dim), vs), vt)
  end
  return BOX0
end

# ///////////////////////////////////////////////////////////
function MAP(fn)
  function MAP0(pol::Hpc)
    if isa(fn, Tuple) || isa(fn, Vector)
      return MapFn(pol, p -> [f(p) for f in fn])
    else
      return MapFn(pol, fn)
    end
  end
  return MAP0
end


# /////////////////////////////////////////////////////////////////
function CIRCLE_POINTS(R, N)
  return [[R * cos(i * 2 * pi / N), R * sin(i * 2 * pi / N)] for i in 0:N-1]
end

# ///////////////////////////////////////////////////////////
function CIRCUMFERENCE(R)
  function CIRCUMFERENCE1(N)
    domain = INTERVALS(2 * pi)(N)
    fn = p -> [R * cos(p[1]), R * sin(p[1])]
    last = domain.childs[1].hulls[end]
    domain.childs[1].hulls[end] = [LEN(domain.childs[1].hulls), 1]
    pop!(domain.childs[1].points)
    return MAP(fn)(domain)
  end
  return CIRCUMFERENCE1
end

# ///////////////////////////////////////////////////////////
function NGON(N)
  return CIRCUMFERENCE(1)(N)
end

# ///////////////////////////////////////////////////////////
function RING(radius::Vector{Float64})
  R1, R2 = radius
  function RING0(subds)
    N, M = subds
    domain = Translate(POWER([INTERVALS(2 * pi)(N), INTERVALS(R2 - R1)(M)]), [0.0, R1])
    fun = p -> [p[2] * cos(p[1]), p[2] * sin(p[1])]
    return MAP(fun)(domain)
  end
  return RING0
end

# ///////////////////////////////////////////////////////////
function TUBE(args::Vector{Float64})
  r1, r2, height = args
  function TUBE0(N)
    return Power(RING([r1, r2])([N, 1]), QUOTE([height]))
  end
  return TUBE0
end

# ///////////////////////////////////////////////////////////
function CIRCLE(R::Number)
  function CIRCLE0(subs)
    N, M = subs
    domain = POWER([INTERVALS(2 * pi)(N), INTERVALS(R)(M)])
    fun = p -> [p[2] * cos(p[1]), p[2] * sin(p[1])]
    return MAP(fun)(domain)
  end
  return CIRCLE0
end

# ///////////////////////////////////////////////////////////
function HOLLOWCYL(rmin=1.0, rmax=2.0, angle=2 * pi, height=2 * rmax)
  function HOLLOWCYL0(shape=[36, 1])
    basis = Hpc(RING(rmin, rmax, angle)(shape))
    return Power(basis, INTERVALS(height)(shape[2]))
  end
  return HOLLOWCYL0
end

# ///////////////////////////////////////////////////////////
function SOLIDCYL(rmax=2.0, angle=2 * pi, height=2 * rmax)
  function SOLIDCYL0(shape=[36, 1])
    return HOLLOWCYL(0.0, rmax, angle, height)(shape)
  end
  return SOLIDCYL0
end

# ///////////////////////////////////////////////////////////
function MY_CYLINDER(args::Vector{Float64})
  R, H = args
  function MY_CYLINDER0(N)
    points = CIRCLE_POINTS(R, N)
    circle = MkPol(points, [collect(1:N)])
    return Power(circle, MkPol([[0], [H]], [[1, 2]]))
  end
  return MY_CYLINDER0
end
CYLINDER = MY_CYLINDER

# /////////////////////////////////////////////////////////////
function CONE(args)
  radius, height = args
  function CONE0(N)
    basis = CIRCLE(radius)([N, 1])
    apex = T(3)(height)(SIMPLEX(0))
    return JOIN([basis, apex])
  end
  return CONE0
end

# /////////////////////////////////////////////////////////////
function TRUNCONE(args)
  R1, R2, H = args
  function TRUNCONE0(N)
    domain = Power(QUOTE([2 * pi / N for i in 1:N]), QUOTE([1]))
    fn = p -> [
      (R1 + p[2] * (R2 - R1)) * cos(p[1]),
      (R1 + p[2] * (R2 - R1)) * sin(p[1]),
      H * p[2]
    ]
    return MAP(fn)(domain)
  end
  return TRUNCONE0
end

# /////////////////////////////////////////////////////////////
function DODECAHEDRON()
  a = 1.0 / (sqrt(3.0))
  g = 0.5 * (sqrt(5.0) - 1)

  top = MKPOL(
    [[1 - g, 1, 0 - g], [1 + g, 1, 0 - g]],
    [[1, 2]],
    [[1]]
  )
  basis = EMBED(1)(CUBOID([2.0, 2.0]))
  roof = T([1, 2, 3])([-1, -1, -1])(JOIN([basis, top]))
  roofpair = STRUCT([roof, R([2, 3])(pi), roof])
  return S([1, 2, 3])([a, a, a])(STRUCT([
    Cube(3, -1, +1),
    roofpair,
    R([1, 3])(pi / 2), R([1, 2])(pi / 2),
    roofpair,
    R([1, 2])(pi / 2), R([2, 3])(pi / 2),
    roofpair
  ]))
end

# /////////////////////////////////////////////////////////////
function ICOSAHEDRON()
  g = 0.5 * (sqrt(5) - 1)
  b = 2.0 / (sqrt(5 * sqrt(5)))
  rectx = T([1, 2])([-g, -1])(CUBOID([2 * g, 2]))
  recty = R([1, 3])(pi / 2)(R([1, 2])(pi / 2)(rectx))
  rectz = R([2, 3])(pi / 2)(R([1, 2])(pi / 2)(rectx))
  return S([1, 2, 3])([b, b, b])(JOIN([rectx, recty, rectz]))
end

function TETRAHEDRON()
  return JOIN([T(3)(-1.0 / 3.0)(NGON(3)), MK([0, 0, 1])])
end

function POLYPOINT(points)
  return MkPol(points, [[i] for i in 1:length(points)])
end

function POLYLINE(points)
  return MkPol(points, [[i, i + 1] for i in 1:length(points)-1])
end

function TRIANGLESTRIPE(points)
  cells = [i % 2 == 0 ? [i, i + 1, i + 2] : [i + 1, i, i + 2] for i in 1:length(points)-2]
  return MkPol(points, cells)
end

function TRIANGLEFAN(points)
  cells = [[1, i - 1, i] for i in 2:length(points)]
  return MkPol(points, cells)
end

function MIRROR(D)
  function MIRROR0(pol)
    return STRUCT([S(D)(-1)(pol), pol])
  end
  return MIRROR0
end


# /////////////////////////////////////////////////////////////
function POLYMARKER(type::Int, MARKERSIZE::Float64=0.1)
  A, B = MARKERSIZE, -MARKERSIZE
  marker0 = MkPol([[A], [0], [0], [A], [B], [0], [0], [B]], [[0 + 1, 1 + 1], [1 + 1, 2 + 1], [2 + 1, 3 + 1], [3 + 1, 0 + 1]])
  marker1 = MkPol([[A], [A], [B], [A], [B], [B], [A], [B]], [[0 + 1, 2 + 1], [1 + 1, 3 + 1]])
  marker2 = MkPol([[A], [A], [B], [A], [B], [B], [A], [B]], [[0 + 1, 1 + 1], [1 + 1, 2 + 1], [2 + 1, 3 + 1], [3 + 1, 0 + 1]])
  marker3 = STRUCT([marker0, marker1])
  marker4 = STRUCT([marker0, marker2])
  marker5 = STRUCT([marker1, marker2])
  marker = [marker0, marker1, marker2, marker3, marker4, marker5][mod(type, 6)+1]
  function POLYMARKER_POINTS(points)
    dim = length(points[1])
    axis = collect(1:dim)
    return Struct([T(axis)(point)(marker) for point in points])
  end
  return POLYMARKER_POINTS
end

# /////////////////////////////////////////////////////////////
function BEZIER(U)
  function BEZIER0(controldata_fn)
    N = length(controldata_fn) - 1
    function map_fn(point)
      t = U(point)
      controldata = [isa(fun, Function) ? fun(point) : fun for fun in controldata_fn]
      ret = [0.0 for i in 1:length(controldata[1])]
      for I in 0:N
        weight = CHOOSE([N, I]) * ((1 - t)^(N - I)) * (t^I)
        for K in 1:length(ret)
          ret[K] += weight * controldata[I+1][K]
        end
      end
      return ret
    end
    return map_fn
  end
  return BEZIER0
end

# /////////////////////////////////////////////////////////////
function BEZIERCURVE(controlpoints)
  return BEZIER(S1)(controlpoints)
end

# /////////////////////////////////////////////////////////////
function COONSPATCH(args)
  su0_fn, su1_fn, s0v_fn, s1v_fn = args
  function map_fn(point)
    u, v = point
    su0 = isa(su0_fn, Function) ? su0_fn(point) : su0_fn
    su1 = isa(su1_fn, Function) ? su1_fn(point) : su1_fn
    s0v = isa(s0v_fn, Function) ? s0v_fn(point) : s0v_fn
    s1v = isa(s1v_fn, Function) ? s1v_fn(point) : s1v_fn
    ret = [0.0 for i in 1:length(su0)]
    for K in 1:length(ret)
      ret[K] = (1 - u) * s0v[K] + u * s1v[K] + (1 - v) * su0[K] + v * su1[K] + (1 - u) * (1 - v) * s0v[K] + (1 - u) * v * s0v[K] + u * (1 - v) * s1v[K] + u * v * s1v[K]
    end
    return ret
  end
  return map_fn
end

# /////////////////////////////////////////////////////////////
function RULEDSURFACE(args)
  alpha_fn, beta_fn = args
  function map_fn(point)
    u, v = point
    alpha, beta = alpha_fn(point), beta_fn(point)
    ret = [0.0 for i in 1:length(alpha)]
    for K in 1:length(ret)
      ret[K] = alpha[K] + v * beta[K]
    end
    return ret
  end
  return map_fn
end

# /////////////////////////////////////////////////////////////
function PROFILEPRODSURFACE(args)
  profile_fn, section_fn = args
  function map_fun(point)
    u, v = point
    profile, section = profile_fn(point), section_fn(point)
    ret = [profile[1] * section[1], profile[1] * section[2], profile[3]]
    return ret
  end
  return map_fun
end

# /////////////////////////////////////////////////////////////
function ROTATIONALSURFACE(args)
  profile = args
  function map_fn(point)
    u, v = point
    f, h, g = profile(point)
    ret = [f * cos(v), f * sin(v), g]
    return ret
  end
  return map_fn
end

# /////////////////////////////////////////////////////////////
function CYLINDRICALSURFACE(args)
  alpha_fun = args[1]
  beta_fun = CONS(AA(K)(args[2]))
  return RULEDSURFACE([alpha_fun, beta_fun])
end

# /////////////////////////////////////////////////////////////
function CONICALSURFACE(args)
  apex = args[1]
  alpha_fn = point -> apex
  beta_fn = point -> [args[2](point)[i] - apex[i] for i in 1:length(apex)]
  return RULEDSURFACE([alpha_fn, beta_fn])
end

# /////////////////////////////////////////////////////////////
function CUBICHERMITE(U)
  function CUBICHERMITE0(args)
    p1_fn, p2_fn, s1_fn, s2_fn = args
    function map_fn(point)
      u = U(point)
      u2 = u * u
      u3 = u2 * u
      p1, p2, s1, s2 = [isa(f, Function) ? f(point) : f for f in [p1_fn, p2_fn, s1_fn, s2_fn]]
      ret = [0.0 for i in 1:length(p1)]
      for i in 1:length(ret)
        ret[i] = (2 * u3 - 3 * u2 + 1) * p1[i] + (-2 * u3 + 3 * u2) * p2[i] + (u3 - 2 * u2 + u) * s1[i] + (u3 - u2) * s2[i]
      end
      return ret
    end
    return map_fn
  end
  return CUBICHERMITE0
end

# /////////////////////////////////////////////////////////////
function HERMITE(args)
  P1, P2, T1, T2 = args
  return CUBICHERMITE(S1)([P1, P2, T1, T2])
end

# /////////////////////////////////////////////////////////////
function Q(H::Number)
  return MkPol([[0], [H]], [[1, 2]])
end

# /////////////////////////////////////////////////////////////
function EXTRUDE(args)
  __N, Pol, H = args
  return Power(Pol, Q(H))
end

# /////////////////////////////////////////////////////////////
function MULTEXTRUDE(P)
  function MULTEXTRUDE0(H)
    return Power(P, Q(H))
  end
  return MULTEXTRUDE0
end

# /////////////////////////////////////////////////////////////
function PROJECT(M)
  function PROJECT0(POL)
    vertices, cells, pols = UKPOL(POL)
    vertices = [vert[1:end-M] for vert in vertices]
    return MKPOL(vertices, cells, pols)
  end
  return PROJECT0
end

# /////////////////////////////////////////////////////////////
function SPLITCELLS(scene)
  vertices, cells, pols = UKPOL(scene)
  ret = []
  for c in cells
    push!(ret, MKPOL(vertices, [c], [[1]]))
  end
  return ret
end

# /////////////////////////////////////////////////////////////
function EXTRACT_WIRES(scene)
  return SPLITCELLS(SKEL_1(scene))
end

SPLITPOLS = SPLITCELLS


# /////////////////////////////////////////////////////////////
function PERMUTAHEDRON(d)
  vertices = ToFloat64(PERMUTATIONS(collect(1:d+1)))
  center = MEANPOINT(vertices)
  cells = [collect(1:length(vertices))]
  object = MKPOL(vertices, cells, [[1]])
  object = Translate(object, [-coord for coord in center])
  for i in 1:d
    object = R([i, d + 1])(pi / 4)(object)
  end
  object = PROJECT(1)(object)
  return object
end

# /////////////////////////////////////////////////////////////
function STAR(N)
  function CIRCLEPOINTS(STARTANGLE)
    function CIRCLEPOINTS1(R)
      function CIRCLEPOINTS0(N)
        return [COMP([CONS([RAISE(PROD)([K(R), COS]), RAISE(PROD)([K(R), SIN])]), RAISE(SUM)([ID, K(STARTANGLE)])]) for STARTANGLE in COMP([COMP([AA(RAISE(PROD)), TRANS]), CONS([K(collect(1:N)), DIESIS(N)])])((2 * pi / N))]
      end
      return CIRCLEPOINTS0
    end
    return CIRCLEPOINTS1
  end
  return COMP([COMP([TRIANGLEFAN, CAT]), TRANS])([CIRCLEPOINTS(0)(1)(N), CIRCLEPOINTS((pi / N))(2.5)(N)])
end


# /////////////////////////////////////////////////////////////
function SCHLEGEL2D(D)
  function map_fn(point)
    return [D * point[1] / point[3], D * point[2] / point[3]]
  end
  return MAP(map_fn)
end

# /////////////////////////////////////////////////////////////
function SCHLEGEL3D(D)
  function map_fn(point)
    return [D * point[1] / point[4], D * point[2] / point[4], D * point[3] / point[4]]
  end
  return MAP(map_fn)
end

# /////////////////////////////////////////////////////////////
function FINITECONE(pol)
  point = [0.0 for i in 1:RN(pol)]
  return JOIN([pol, MK(point)])
end

# /////////////////////////////////////////////////////////////
function PRISM(HEIGHT)
  function PRISM0(BASIS)
    return Power(BASIS, QUOTE([HEIGHT]))
  end
  return PRISM0
end

# /////////////////////////////////////////////////////////////
function CROSSPOLYTOPE(D)
  points = []
  for i in 1:D
    point_pos = [0 for x in 1:D]
    point_pos[i] = +1
    point_neg = [0 for x in 1:D]
    point_neg[i] = -1
    push!(points, point_pos, point_neg)
  end
  cells = [collect(1:D*2)]
  pols = [[1]]
  return MKPOL(points, cells, pols)
end

# /////////////////////////////////////////////////////////////
function OCTAHEDRON()
  return CROSSPOLYTOPE(2)
end

# /////////////////////////////////////////////////////////////
function ROTN(args)
  alpha, N = args
  N = UNITVECT(N)
  QX = UNITVECT(VECTPROD([[0, 0, 1], N]))
  QZ = UNITVECT(N)
  QY = VECTPROD([QZ, QX])
  Q = MATHOM([QX, QY, QZ])
  ISUP = COMP([AND, CONS([COMP([C(EQ)(0), S1]), COMP([C(EQ)(0), S2]), COMP([COMP([NOT, C(EQ)(0)]), S3])])])
  if N[1] == 0 && N[2] == 0
    return R([1, 2])(alpha)
  else
    return COMP([MAT(TRANS(Q)), R([1, 2])(alpha), MAT(Q)])
  end
end

# /////////////////////////////////////////////////////////////
function MKVERSORK()
  return TOP([CYLINDER([1.0 / 100.0, 7.0 / 8.0])(6), CONE([1.0 / 16.0, 1.0 / 8])(8)])
end

function MKVECTOR(P1)
  function MKVECTOR0(P2)
    TR = T([1, 2, 3])(P1)
    U = VECTDIFF([P2, P1])
    ALPHA = acos(INNERPROD([[0, 0, 1], UNITVECT(U)]))
    B = VECTNORM(U)
    SC = S([1, 2, 3])([B, B, B])
    N = VECTPROD([[0, 0, 1], U])
    ROT = ROTN([ALPHA, N])
    return COMP([COMP([TR, ROT]), SC])(MKVERSORK())
  end
  return MKVECTOR0
end



# /////////////////////////////////////////////////////////////
function CUBICUBSPLINE(domain)
  function CUBICUBSPLINE0(args)
    q1_fn, q2_fn, q3_fn, q4_fn = args
    function map_fn(point)
      u = S1(point)
      u2 = u * u
      u3 = u2 * u
      q1, q2, q3, q4 = [isa(f, Function) ? f(point) : f for f in [q1_fn, q2_fn, q3_fn, q4_fn]]
      ret = [0.0 for x in 1:length(q1)]
      for i in 1:length(ret)
        ret[i] = (1.0 / 6.0) * ((-u3 + 3 * u2 - 3 * u + 1) * q1[i] + (3 * u3 - 6 * u2 + 4) * q2[i] + (-3 * u3 + 3 * u2 + 3 * u + 1) * q3[i] + (u3) * q4[i])
      end
      return ret
    end
    return MAP(map_fn)(domain)
  end
  return CUBICUBSPLINE0
end

# //////////////////////////////////////////////////////////////
function CUBICCARDINAL(domain, h=1)
  function CUBICCARDINAL0(args)
    q1_fn, q2_fn, q3_fn, q4_fn = args
    function map_fn(point)
      u = S1(point)
      u2 = u * u
      u3 = u2 * u
      q1, q2, q3, q4 = [isa(f, Function) ? f(point) : f for f in [q1_fn, q2_fn, q3_fn, q4_fn]]
      ret = [0.0 for i in 1:length(q1)]
      for i in 1:length(ret)
        ret[i] = (-h * u3 + 2 * h * u2 - h * u) * q1[i] + ((2 - h) * u3 + (h - 3) * u2 + 1) * q2[i] + ((h - 2) * u3 + (3 - 2 * h) * u2 + h * u) * q3[i] + (h * u3 - h * u2) * q4[i]
      end
      return ret
    end
    return MAP(map_fn)(domain)
  end
  return CUBICCARDINAL0
end

# //////////////////////////////////////////////////////////////
function SPLINE(curve)
  function SPLINE0(points)
    ret = Vector{Hpc}()
    for i in 1:length(points)-4+1
      P = points[i:i+4-1]
      push!(ret, curve(P))
    end
    return Struct(ret)
  end
  return SPLINE0
end

# //////////////////////////////////////////////////////////////
function JOINTS(curve)
  knotzero = MK([0.0])
  function JOINTS0(points)
    points, cells, pols = UKPOL(SPLINE(curve(knotzero)))
    return POLYMARKER(2)(points)
  end
  return JOINTS0
end

# //////////////////////////////////////////////////////////////
function BERNSTEINBASIS(U)
  function BERNSTEIN0(N)
    function BERNSTEIN1(I)
      function map_fn(point)
        t = U(point)
        ret = CHOOSE([N, I]) * ((1 - t)^(N - I)) * (t^I)
        return ret
      end
      return map_fn
    end
    return [BERNSTEIN1(I) for I in 0:N]
  end
  return BERNSTEIN0
end

# //////////////////////////////////////////////////////////////
function TENSORPRODSURFACE(args)
  ubasis, vbasis = args
  function TENSORPRODSURFACE0(controlpoints_fn)
    controlpoints_fn = ToFloat64(controlpoints_fn)
    function map_fn(point)
      u, v = point
      U = [f([u]) for f in ubasis]
      V = [f([v]) for f in vbasis]
      controlpoints = [isa(f, Function) ? f(point) : f for f in controlpoints_fn]
      target_dim = length(controlpoints[1][1])
      ret = [0.0 for x in 1:target_dim]
      for i in 1:length(ubasis)
        for j in 1:length(vbasis)
          for M in 1:target_dim
            ret[M] += U[i] * V[j] * controlpoints[i][j][M]
          end
        end
      end
      return ret
    end
    return map_fn
  end
  return TENSORPRODSURFACE0
end

# //////////////////////////////////////////////////////////////
function BILINEARSURFACE(controlpoints)
  return TENSORPRODSURFACE([BERNSTEINBASIS(S1)(1), BERNSTEINBASIS(S1)(1)])(controlpoints)
end

# //////////////////////////////////////////////////////////////
function BIQUADRATICSURFACE(controlpoints)

  function u0(point)
    u = S1(point)
    return 2 * u * u - u
  end
  function u1(point)
    u = S1(point)
    return 4 * u - 4 * u * u
  end
  function u2(point)
    u = S1(point)
    return 2 * u * u - 3 * u + 1
  end
  basis = [u0, u1, u2]
  return TENSORPRODSURFACE([basis, basis])(controlpoints)
end

# //////////////////////////////////////////////////////////////
function HERMITESURFACE(controlpoints)
  function H0(point)
    u = S1(point)
    u2 = u * u
    u3 = u2 * u
    return u3 - u2
  end
  function H1(point)
    u = S1(point)
    u2 = u * u
    u3 = u2 * u
    return u3 - 2 * u2 + u
  end
  function H2(point)
    u = S1(point)
    u2 = u * u
    u3 = u2 * u
    return 3 * u2 - 2 * u3
  end
  function H3(point)
    u = S1(point)
    u2 = u * u
    u3 = u2 * u
    return 2 * u3 - 3 * u2 + 1
  end

  basis = [H3, H2, H1, H0]
  return TENSORPRODSURFACE([basis, basis])(controlpoints)
end

# //////////////////////////////////////////////////////////////
function BEZIERSURFACE(controlpoints)
  M = length(controlpoints) - 1
  N = length(controlpoints[1]) - 1
  return TENSORPRODSURFACE([BERNSTEINBASIS(S1)(M), BERNSTEINBASIS(S1)(N)])(controlpoints)
end

# //////////////////////////////////////////////////////////////
function TENSORPRODSOLID(args)
  ubasis, vbasis, wbasis = args
  function TENSORPRODSOLID0(controlpoints_fn)
    controlpoints_fn = ToFloat64(controlpoints_fn)
    function map_fn(point)
      u, v, w = point
      U = [f([u]) for f in ubasis]
      V = [f([v]) for f in vbasis]
      W = [f([w]) for f in wbasis]
      controlpoints = [isa(f, Function) ? f(point) : f for f in controlpoints_fn]
      target_dim = length(controlpoints[1][1][1])
      ret = [0.0 for x in 1:target_dim]
      for i in 1:length(ubasis)
        for j in 1:length(vbasis)
          for k in 1:length(wbasis)
            for M in 1:target_dim
              ret[M] += U[i] * V[j] * W[k] * controlpoints[M][i][j][k]
            end
          end
        end
      end
      return ret
    end
    return map_fn
  end
  return TENSORPRODSOLID0
end

# //////////////////////////////////////////////////////////////
function BEZIERMANIFOLD(degrees)
  basis = [BERNSTEINBASIS(S1)(d) for d in degrees]
  return TENSORPRODSOLID(basis)
end

function LOCATE(args)
  pol, a, distances = args
  ret = []
  for d in distances
    push!(ret, T(a)(d), pol)
  end
  return STRUCT(ret)
end

function RIF(size)
  thin = 0.01 * size
  x = COLOR(RED)(CUBOID([size, thin, thin]))
  y = COLOR(GREEN)(CUBOID([thin, size, thin]))
  z = COLOR(BLUE)(CUBOID([thin, thin, size]))
  return STRUCT([x, y, z])
end

# //////////////////////////////////////////////////////////////
function FRACTALSIMPLEX(D)
  function FRACTALSIMPLEX0(N)
    mkpols = COMP([COMP([COMP([COMP([STRUCT, AA(MKPOL)]), AA(AL)]), DISTR]), CONS([ID, K([[FROMTO([1, D + 1])], [[1]]])])])
    function COMPONENT(args)
      i, seq = args
      firstseq = seq[1:i-1]
      pivot = seq[i-1]
      lastseq = seq[i:length(seq)]
      firstpart = AA(MEANPOINT)(DISTR([firstseq, pivot]))
      lastpart = AA(MEANPOINT)(DISTR([lastseq, pivot]))
      return CAT([firstpart, [pivot], lastpart])
    end
    expand = COMP([COMP([AA(COMPONENT), DISTR]), CONS([COMP([INTSTO, LEN]), ID])])
    splitting = COMP([COMP, DIESIS(N)])(COMP([CAT, AA(expand)]))
    return COMP([COMP([COMP([COMP([mkpols, splitting]), CONS([S1])])])])(UKPOL(SIMPLEX(D)))
  end
  return FRACTALSIMPLEX0
end


function MESH(seq)
  return INSL(RAISE(PROD))([QUOTE(i) for i in seq])
end

function NU_GRID(data)
  polylines = [POLYLINE(i) for i in data]
  return INSL(RAISE(PROD))(polylines)
end


function SEGMENT(sx)
  function SEGMENT0(args)
    N = length(args[1])
    A, B = args
    P0 = A
    P1 = [A[i] + (B[i] - A[i]) * sx for i in 1:N]
    return POLYLINE([P0, P1])
  end
  return SEGMENT0
end

# //////////////////////////////////////////////////////////////
function SOLIDIFY(pol)
  mybox = box(pol)
  min = mybox.p1[1]
  max = mybox.p2[1]
  siz = max - min
  far_point = max + siz * 100
  function InftyProject(pol)
    verts, cells, pols = UKPOL(pol)
    verts = [[far_point] + v[2:end] for v in verts]
    return MKPOL(verts, cells, pols)
  end
  function IsFull(pol)
    return DIM(pol) == RN(pol)
  end
  ret = SPLITCELLS(pol)
  ret = [JOIN([pol, InftyProject(pol)]) for pol in ret]
  objs = convert(Vector{Hpc}, FILTER(IsFull)(ret))
  return XOR(objs)
end

# //////////////////////////////////////////////////////////////
function EXTRUSION(angle)
  function EXTRUSION1(height)
    function EXTRUSION0(pol)
      dim = DIM(pol)
      cells = SPLITCELLS(SKELETON(dim)(pol))
      slice = [EMBED(1)(c) for c in cells]
      tensor = COMP([T(dim + 1)(1.0 / height), R([dim - 1, dim])(angle / height)])
      layer = STRUCT([JOIN([p, tensor(p)]) for p in slice])
      return COMP([COMP([STRUCT, CAT]), DIESIS(height)])([layer, tensor])
    end
    return EXTRUSION0
  end
  return EXTRUSION1
end

# //////////////////////////////////////////////////////////////
function EX(args)
  x1, x2 = args
  function EX0(pol)
    dim = DIM(pol)
    return T(dim + 1)(x1)(S(dim + 1)(x2 - x1)(EXTRUSION(0.0)(1.0)(pol)))
  end
  return EX0
end

# //////////////////////////////////////////////////////////////
function LEX(args)
  x1, x2 = args
  function LEX0(pol)
    function SHEARTENSOR(A)
      function SHEARTENSOR0(POL)
        dim = DIM(POL)
        newrow = K((AR([CAT([[0, 1], DIESIS((dim - 2))(0)]), A])))
        update = (COMP([CONS, CAT]))([[S1, newrow], AA(SEL)((FROMTO([3, dim + 1])))])
        matrix = update(IDNT(dim + 1))
        return (MAT(matrix))(POL)
      end
      return SHEARTENSOR0
    end
    ret = EXTRUSION(0)(1)(pol)
    ret = SHEARTENSOR(x2 - x1)(ret)
    ret = S(DIM(pol) + 1)(x2 - x1)(ret)
    ret = T(DIM(pol) + 1)(x1)(ret)
    return ret
  end
  return LEX0
end

# //////////////////////////////////////////////////////////////
function SEX(args)
  x1, x2 = args
  function SEX1(height)
    function SEX0(pol)
      dim = DIM(pol)
      ret = EXTRUSION(x2 - x1)(height)(pol)
      ret = S(dim + 1)(x2 - x1)(ret)
      ret = R([dim, dim - 1])(x1)(ret)
      return ret
    end
    return SEX0
  end
  return SEX1
end

# //////////////////////////////////////////////////////////////
function UKPOLF(pol)
  error("not implemented")
end

# //////////////////////////////////////////////////////////////
function POLAR(pol, precision=1e-6)
  faces, cells, pols = UKPOLF(pol)
  for i in 1:length(faces)
    mod = -1 * faces[i][1]
    if abs(mod) < precision
      mod = 1
    end
    faces[i] = [value / mod for value in faces[i][2:end]]
  end
  return MKPOL(faces, cells, pols)
end

# //////////////////////////////////////////////////////////////
function SWEEP(v)
  function SWEEP0(pol)
    ret = Power(pol, QUOTE([1]))
    mat = IDNT(length(v) + 2)
    for i in 1:length(v)
      mat[i+1, length(v)+1] = v[i]
    end
    ret = MAT(mat)(ret)
    return PROJECT(1)(ret)
  end
  return SWEEP0
end

# //////////////////////////////////////////////////////////////
function MINKOWSKI(vects)
  function MINKOWSKI0(pol)
    ret = pol
    for I in length(vects):-1:1
      ret = SWEEP(vects[I])(ret)
    end
    return ret
  end
  return MINKOWSKI0
end

# //////////////////////////////////////////////////////////////
function OFFSET(v)
  function OFFSET0(pol)
    ret = pol
    for I in 1:length(v)
      shear = [[J != I ? 0.0 : v[I] for J in 1:length(v)]; [0.0 for J in 2:I]]
      mat = IDNT(length(shear) + 2)
      for J in 1:length(shear)
        mat[J+1, length(shear)+2] = shear[J]
      end
      ret = MAT(mat)((Power(ret, QUOTE([1]))))
    end
    return PROJECT(length(v))(ret)
  end
  return OFFSET0
end

# //////////////////////////////////////////////////////////////
function THINSOLID(surface, delta=1e-4)
  function map_fn(point)
    u, v, w = point
    P0 = surface([u, v])
    PX = surface([u + delta, v])
    PY = surface([u, v + delta])
    GX = [PX[i] - P0[i] for i in 1:3]
    GY = [PY[i] - P0[i] for i in 1:3]
    normal = UNITVECT(VECTPROD([GX, GY]))
    ret = [P0[i] + w * normal[i] for i in 1:3]
    return ret
  end
  return map_fn
end

# //////////////////////////////////////////////////////////////
function PLANE(args)
  p0, p1, p2 = args
  v1 = VECTDIFF([p1, p0])
  v2 = VECTDIFF([p2, p0])
  side1 = VECTNORM(v1)
  side2 = VECTNORM(v2)
  normal = UNITVECT(VECTPROD([v1, v2]))
  axis = VECTPROD([[0, 0, 1], normal])
  angle = acos(INNERPROD([[0, 0, 1], normal]))
  geometry = T([1, 2, 3])(p0)(ROTN([angle, axis])(T([1, 2])([-1 * side1, -1 * side2])(CUBOID([2 * side1, 2 * side2]))))
  return [normal, p0, geometry]
end

# //////////////////////////////////////////////////////////////
function RATIONALBEZIER(controlpoints_fn)
  controlpoints_fn = ToFloat64(controlpoints_fn)
  degree = length(controlpoints_fn) - 1
  basis = BERNSTEINBASIS(S1)(degree)
  function map_fn(point)
    controlpoints = [isa(f, Function) ? f(point) : f for f in controlpoints_fn]
    target_dim = length(controlpoints[1])
    ret = [0.0 for i in 1:target_dim]
    for i in 1:length(basis)
      coeff = basis[i](point)
      for M in 1:target_dim
        ret[M] += coeff * controlpoints[i][M]
      end
    end
    last = ret[end]
    if last != 0
      ret = [value / last for value in ret]
    end
    ret = ret[1:end-1]
    return ret
  end
  return map_fn
end

# //////////////////////////////////////////////////////////////
function ELLIPSE(args::Vector{Float64})
  A, B = args
  function ELLIPSE0(N::Int)
    C = 0.5 * sqrt(2)
    mapping = RATIONALBEZIER([[A, 0.0, 1.0], [A * C, B * C, C], [0.0, B, 1.0]])
    quarter = MAP(mapping)((INTERVALS(1.0)(N)))
    half = STRUCT([quarter, S(2)(-1.0)(quarter)])
    return STRUCT([half, S(1)(-1.0)(half)])
  end
  return ELLIPSE0
end

# //////////////////////////////////////////////////////////////
function CURVE_NORMAL(curve)
  function map_fn(point)
    xu, yu = curve(point)
    mod2 = xu * xu + yu * yu
    den = mod2 > 0 ? sqrt(mod2) : 0
    return [-yu / den, xu / den]
  end
  return map_fn
end

# //////////////////////////////////////////////////////////////
function DERBEZIER(controlpoints_fn)
  controlpoints_fn = ToFloat64(controlpoints_fn)
  degree = length(controlpoints_fn) - 1

  function DERBERNSTEIN(N)
    function DERBERNSTEIN0(I)
      function map_fn(point)
        t = S1(point)
        return CHOOSE([N, I]) * t^(I - 1) * (1 - t)^(N - I - 1) * (I - N * t)
      end
      return map_fn
    end
    return DERBERNSTEIN0
  end

  basis = [DERBERNSTEIN(degree)(i) for i in 1:degree+1]
  function map_fn(point)
    controlpoints = [isa(f, Function) ? f(point) : f for f in controlpoints_fn]
    target_dim = length(controlpoints[1])
    ret = [0.0 for i in 1:target_dim]
    for i in 1:length(basis)
      coeff = basis[i](point)
      for M in 1:target_dim
        ret[M] += coeff * controlpoints[i][M]
      end
    end
    return ret
  end
  return map_fn
end

# //////////////////////////////////////////////////////////////
function BEZIERSTRIPE(args)
  controlpoints, width, n = args
  bezier = BEZIERCURVE(controlpoints)
  normal = CURVE_NORMAL(DERBEZIER(controlpoints))
  function map_fn(point)
    u, v = point
    bx, by = bezier(point)
    nx, ny = normal(point)
    ret = [bx + v * nx, by + v * ny]
    return ret
  end
  domain = S(2)(width)(T(1)(0.00001)(Power(INTERVALS(1.0)(n), INTERVALS(1.0)(1))))
  return MAP(map_fn)(domain)
end

# //////////////////////////////////////////////////////////////
function BSPLINE(degree::Int)
  function BSPLINE0(knots)
    function BSPLINE1(points_fn)
      n = length(points_fn) - 1
      m = length(knots) - 1
      k = degree + 1
      T = knots
      tmin, tmax = T[k], T[n+1]
      @assert length(knots) == (n + k + 1)

      function N(i, k, t)

        if k == 1
          if (t >= T[i] && t < T[i+1]) || (t == tmax && t >= T[i] && t <= T[i+1])
            return 1
          else
            return 0
          end
        end

        ret = 0
        num1, div1 = t - T[i], T[i+k-1] - T[i]
        if div1 != 0
          ret += (num1 / div1) * N(i, k - 1, t)
        end
        num2, div2 = T[i+k] - t, T[i+k] - T[i+1]
        if div2 != 0
          ret += (num2 / div2) * N(i + 1, k - 1, t)
        end
        return ret
      end

      function map_fn(point)
        t = point[1]

        points = [isa(f, Function) ? f(point) : f for f in points_fn]
        target_dim = length(points[1])
        ret = [0.0 for i in 1:target_dim]
        for i in 1:n+1
          coeff = N(i, k, t)
          for M in 1:target_dim
            ret[M] += points[i][M] * coeff
          end
        end
        return ret
      end
      return map_fn
    end
    return BSPLINE1
  end
  return BSPLINE0
end

# //////////////////////////////////////////////////////////////
function NUBSPLINE(degree, totpoints=80)
  function NUBSPLINE1(knots)
    function NUBSPLINE2(points)
      m = length(knots)
      tmin = min(knots...)
      tmax = max(knots...)
      tsiz = tmax - tmin
      v = [tsiz / float(totpoints - 1) for i in 1:totpoints-1]
      @assert length(v) + 1 == totpoints
      v = [-tmin; v]
      domain = QUOTE(v)
      return MAP(BSPLINE(degree)(knots)(points))(domain)
    end
    return NUBSPLINE2
  end
  return NUBSPLINE1
end

# //////////////////////////////////////////////////////////////
function DISPLAYNUBSPLINE(degree::Int, knots::Vector, points, marker_size::Float64=0.1)
  polymarker_type = 2
  spline = NUBSPLINE(degree, length(knots))(knots)(points)
  spline_view_knots = POLYMARKER(polymarker_type, marker_size)(UKPOL(spline)[1])
  return STRUCT([
    degree > 0 ? NUBSPLINE(degree)(knots)(points) : POLYMARKER(3, marker_size)(points), spline_view_knots, POLYLINE(points), POLYMARKER(1, marker_size)(points)
  ])
end

# //////////////////////////////////////////////////////////////
function RATIONALBSPLINE(degree)
  function RATIONALBSPLINE0(knots)
    function RATIONALBSPLINE1(points)
      bspline = BSPLINE(degree)(knots)(points)
      function map_fn(point)
        ret = bspline(point)
        last = ret[end]
        if last != 0
          ret = [value / last for value in ret]
        end
        ret = ret[1:end-1]
        return ret
      end
      return map_fn
    end
    return RATIONALBSPLINE1
  end
  return RATIONALBSPLINE0
end

# //////////////////////////////////////////////////////////////
function NURBSPLINE(degree, totpoints=80)
  function NURBSPLINE1(knots)
    function NURBSPLINE2(points)
      m = length(knots)
      tmin = min(knots...)
      tmax = max(knots...)
      tsiz = tmax - tmin
      v = [tsiz / float(totpoints - 1) for i in 1:totpoints-1]
      @assert length(v) + 1 == totpoints
      v = [-tmin; v]
      domain = QUOTE(v)
      return MAP(RATIONALBSPLINE(degree)(knots)(points))(domain)
    end
    return NURBSPLINE2
  end
  return NURBSPLINE1
end

# //////////////////////////////////////////////////////////////
function DISPLAYNURBSPLINE(args, marker_size=0.1)
  degree, knots, points = args
  spline_view_knots = POLYMARKER(2, marker_size)(UKPOL(NURBSPLINE(degree, length(knots))(knots)(points))[1])
  return STRUCT([
    degree > 0 ?
    NURBSPLINE(degree)(knots)(points) :
    POLYMARKER(3, marker_size)(points), spline_view_knots, POLYLINE(points), POLYMARKER(1, marker_size)(points)
  ])
end


# //////////////////////////////////////////////////////////////
function COLOR(C)
  function COLOR0(hpc::Hpc)
    return PROPERTIES(hpc, Properties("face_color" => Point4d(C[1], C[2], C[3], length(C) >= 4 ? C[4] : 1.0)))
  end
  return COLOR0
end

# //////////////////////////////////////////////////////////////
function PROPERTIES(hpc::Hpc, properties::Properties)
  ret = STRUCT([hpc])
  ret.properties = copy(hpc.properties)
  for (key, value) in properties
    ret.properties[key] = value
  end
  return ret
end


# //////////////////////////////////////////////////////////////
function ICOSPHERE(obj::Hpc=ICOSAHEDRON())::Hpc  # obj = Hpc
  W = LAR(obj).V
  EV = LAR(obj).C[:EV]
  W = [W[:, k] for k = 1:size(W, 2)]
  V = [(W[v1] + W[v2]) ./ 2 for (v1, v2) in EV]
  r1 = sqrt(sum(W[1] .^ 2))
  s1 = sqrt(sum(V[1] .^ 2))
  CONVEXHULL([W; V * (r1 / s1)])
end

# /////////////////////////////////////////////////////////////
function icosphere(level::Int64=1)
  obj0 = ICOSAHEDRON()
  if level == 0
    return obj0
  end
  if level == 1
    return obj1 = ICOSPHERE(obj0)
  end
  if level >= 2
    obj1 = ICOSPHERE(obj0)
    return ICOSPHERE(obj1)
  end
end


# //////////////////////////////////////////////////////////////////////////////
function SPHERE(radius=1.0::Number)
  function SPHERE0(subds=[16, 32]::Vector{Int})
    N, M = subds
    domain = T(1, 2)(-pi / 2, -pi)(Power(INTERVALS(pi)(N), INTERVALS(2 * pi)(M)))
    fx = p -> radius * (-cos(p[1])) * sin(p[2])
    fy = p -> radius * cos(p[1]) * cos(p[2])
    fz = p -> radius * sin(p[1])
    return MAP([fx, fy, fz])(domain)
  end
  return SPHERE0
end

# //////////////////////////////////////////////////////////////////////////////
function TORUS(radii=[1.0, 2]::Vector)
  r1, r2 = radii
  function TORUS0(subds=[16, 32]::Vector{Int})
    N, M = subds
    a = 0.5 * (r2 - r1)
    c = 0.5 * (r1 + r2)
    domain = Power(INTERVALS(2 * pi)(N), INTERVALS(2 * pi)(M))
    fx = p -> (c + a * cos(p[2])) * cos(p[1])
    fy = p -> (c + a * cos(p[2])) * sin(p[1])
    fz = p -> a * sin(p[2])
    return MAP([fx, fy, fz])(domain)
  end
  return TORUS0
end

# //////////////////////////////////////////////////////////////////////////////
"""
    GRID1(n::Int)::Hpc
Generate a 1D object of `Hpc` type with `n` unit segments.
```
"""
GRID1(n) = QUOTE(DIESIS(n)(1.0))
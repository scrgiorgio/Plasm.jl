using Plasm

# ////////////////////////////////////////////////////////
function RandomLine(size_min::Float64,size_max::Float64)
  size = size_min+rand()*(size_max-size_min)
  return STRUCT(
    T(1,2)(rand(2)...), 
    S([1,2])([size,size]), 
    R([1,2])(2*pi*rand()),
    Plasm.SQUARE(1)
  )
end

# ////////////////////////////////////////////////////////
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

# make sure I will get consistent random numbers (important for debugging)
import Random
Random.seed!(0)

hpc = STRUCT([RandomLine(2.0,3.0) for I in 1:6])
# VIEW(hpc)

lar = LAR(hpc)
V, EV  = lar.V, lar.C[:EV]
V,FVs,EVs = Plasm.arrange2D(V,EV)

ViewColored(V, EVs, 1.0)
ViewColored(V, EVs, 1.2)

ViewColored(V, FVs, 1.0)
ViewColored(V, FVs, 1.2)








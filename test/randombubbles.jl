using Plasm


# //////////////////////////////////////////////////////////
function RandomBubble()
  vs = rand()
  vt = rand(2)
  return STRUCT(
    T(1,2)(vt...),
    S([1,2])([0.25*vs, 0.25*vs]), 
    CIRCLE(1)([8,1])
  )
end

VIEW(RandomBubble())

# ////////////////////////////////////////////////////////
function ViewColored(V,cells, scale=1.2, line_width=3)
  exploded = explodecells(V, cells, sx=scale, sy=scale, sz=scale)
  v=[]
  for k in eachindex(exploded)
    c = Point4d(Plasm.COLORS[(k-1)%12+1] - (rand(Float64,4)*0.1))
    c[4] = 1.0
    push!(v,PROPERTIES(exploded[k], Dict(
      "line_color" => c, 
      "face_color" => c,
      "line_width" => line_width)))
  end
  VIEW(STRUCT(v))
end

# make sure I will get consistent random numbers (important for debugging)
import Random
Random.seed!(0)

hpc = STRUCT([RandomBubble() for I in 1:50])
VIEW(hpc)

lar = LAR(hpc)
V, EV  = lar.V, lar.C[:EV]

V,FVs,EVs = Plasm.arrange2D(V,EV)

ViewColored(V, EVs, 1.0)
ViewColored(V, EVs, 1.2)
ViewColored(V, FVs, 1.0)
ViewColored(V, FVs, 1.2)








using Plasm
using PlyIO

ply=load_ply("resources/bunny-closed.ply")

x,y,z=[collect(ply["vertex"][it]) for it in ("x","y","z")]
V=ToPoints([Vector{Float32}([X,Y,Z]) for (X,Y,Z) in zip(x,y,z)])

FV=Vector{Vector{Int}}()
for (i,j,k) in ply["face"]["vertex_indices"]
  push!(FV,Vector{Int}([i+1,j+1,k+1]))
end

P=[V,FV]
println(VOLUME(P))
println(SURFACE(P))



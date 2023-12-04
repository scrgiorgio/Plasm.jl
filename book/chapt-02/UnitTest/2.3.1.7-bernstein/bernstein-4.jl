using Plasm
m = 30; EV = [[k, k+1] for k=1:m]
for n=1:4
    meshes = [ GL.GLFrame2 ]
    for k=1:n+1
        V = points2D(Bernstein(n)[k],m) 
        push!(meshes, GL.GLLines(V,EV))
    end
    VIEW(meshes);
end

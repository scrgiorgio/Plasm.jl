function hpcprod( hpc1, hpc2 )
   key1 = findmax(keys, hpc1.C)[2]
   key2 = findmax(keys, hpc2.C)[2]
   (V, cells1) = hpc1.V, values(hpc1.C[key1])
   (W, cells2) = hpc2.V, values(hpc2.C[key2])
   vertices = vertprod(V, W)
   thecells = cellprod(cells1, cells2, V, W)
   cells = [[vertices[v] for v in cell] for cell in thecells]
   verts = hcat(keys(vertices)...); d = size(verts, 1)
   return Lar(verts, Dict(Symbol("c$(d)v") => cells))
end

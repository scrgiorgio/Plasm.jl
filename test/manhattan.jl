using Plasm
import Base.*, Base.splat
*(a::Hpc, b::Hpc) = Power(a, b)

function Manhattan2D()
   global verts = [[0.,0],[3,0],[5,0],[7,0],[8,0],[9.5,1],[10,1.5],[0,3],[3,3],
   [5,3],[7,3],[8,3],[9.5,3],[0,4],[3,4],[5,4],[9.5,4],[12,4],[9.5,5],[10,5],
   [12,5],[0,6],[3,6],[5,6],[0,7],[3,7],[5,7],[9.5,7],[12,7],[9.5,8],
   [12,8],[0,9],[3,9],[5,9],[8,9],[9,9],[12,9],[0,10],[3,10],[5,10],
   [8,10],[9,10],[9.5,10],[10,10],[12,10],[6,11],[7,11],[0,12],[3,12],[9,12],
   [9.5,12],[0,13],[3,13],[6,13],[7,13],[9,13],[9.5,13],[0,14],[3,14],[5,14],
   [8,14],[9,14],[9.5,14],[10,14],[12,14],[0,15],[3,15],[5,15],[8,15],[0,16],
   [6,16],[7,16],[9,17],[9.5,17],[10,17],[12,17],[6,18],[7,18],[9,18],[9.5,18],
   [10,18],[12,18],[2,19],[3,19],[5,19],[8,19],[9,19],[9.5,19],[10,19],[12,19],
   [5,20],[12,20],[7,22],[10,22],[9,6],[12,6],[9,15],[9.5,15],[10,15],[12,15]];
   global cells =[[1,2,9,8],[3,4,11,10],[5,6,13,12],[14,15,23,22],[16,17,19,24],
   [7,18,21,20],[25,26,33,32],[27,95,28,35,34],[95,96,29,28],[30,31,37,36],
   [38,39,49,48],[40,41,47,46],[41,61,55,47],[55,61,60,54],[54,60,40,46],
   [42,43,51,50],[44,45,65,64],[52,53,59,58],[56,57,63,62],[66,67,84,83,70],
   [68,69,72,71],[69,86,78,72],[78,86,85,77],[71,77,85,68],[97,98,74,73],
   [99,100,76,75],[79,80,88,87],[81,82,90,89], [91,92,94,93]];
   model = MKPOL(verts,cells)
   VIEW( model, "Manhattan2D" )
end

function Manhattan3D()
   ManhattanH = [1,3,1,11,1,2,1,1,1,8,15,1,1,1,1,8,1,15,8, 1,2,2,2,2,5,9,1,1,1].*3
   # 29-element Vector{Int64}:
   storeys = CONS(AA(DIESIS)(ManhattanH))(.5)
   # 29-element Vector{Vector{Float64}}:
   pols1D = AA(QUOTE)(storeys)
   # 29-element Vector{Hpc}:
   pols2D = [MKPOL(verts,[cell]) for cell in cells]
   # 29-element Vector{Hpc}:
   pols3D = AA(splat(*))(TRANS([pols2D, pols1D]))
   # 29-element Vector{Hpc}:
   VIEW(STRUCT(pols3D), "Manhattan3D")
end

Manhattan2D()
Manhattan3D()

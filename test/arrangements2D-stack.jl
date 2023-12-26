# ///////////////////////////////////////////////////////
using Plasm 
include("../src/arrangement2D.jl")

function explodecells(V,FVs; sx=1.2,sy=1.2,sz=1.2) 
   ret = []
   for FV in FVs  
   # FV = single cell made by edges or trias or quads
      # remove repeated vert indices
      vertsidx = !(FV==Any[]) ? sort(union(FV...)) : nothing  
      # geometry matrix of cell
      vcell = V[:,vertsidx] 
      # key=vertsidx => value=ordinal n.
      vdict = Dict(zip(vertsidx,1:length(vertsidx))) 
      # barycenter: sum of points in FV cell/their number
      centre = sum(vcell,dims=2)/size(vcell,2)  
      scaled_center = size(centre,1)==2 ? centre .* [sx;sy] : centre .* [sx;sy;sz]
      translation_vector = scaled_center - centre
      # traslated geometry matrix of cell
      cellverts = vcell .+ translation_vector 
      # points of Hpc value
      points = [cellverts[:,k] for k=1:size(cellverts,2)]
      # local topology of chain primitive (array of triangles)
      chain = [[vdict[v] for v in cell] for cell in FV] 
      # single item is MKPOL argument
      push!(ret, MKPOL(points, chain))
   end
   return ret
end

function TestExplode()
   obj=ToLAR(STRUCT(
      Plasm.SQUARE(1), 
      T(1,2)(0.5,0.25), 
      Plasm.SQUARE(1)
   ))

   V = hcat(obj.childs[1].points...)
   EV = obj.childs[1].edges
   V,FVs,EVs = Plasm.arrange2D(V,EV)
  
   begin
      exploded = explodecells(V, FVs, sx=1.2, sy=1.2, sz=1.2)
      v=[]
      for k in 1:length(exploded)
        face_color = Point4d(Plasm.COLORS[(k-1)%12+1] - (rand(Float64,4)*0.1))
        face_color[4] = 1.0
        push!(v,PROPERTIES(exploded[k], Dict(
        "face_color" => face_color, 
        #"line_color" => GL.BLACK, 
        "line_width" => 3)))
      end
      VIEW(STRUCT(v))
    end
  
   begin
      exploded=explodecells(V, EVs, sx=1.2, sy=1.2, sz=1.2)
      v=[]
      for k in 1:length(exploded)
         line_color=Point4d(Plasm.COLORS[(k-1)%12+1] - (rand(Float64,4)*0.1))
         line_color[4]=1.0    
         push!(v,PROPERTIES(exploded[k], Dict("line_color" => line_color, "line_width" => 3)))
      end
      VIEW(STRUCT(v))
   end
end

TestExplode()



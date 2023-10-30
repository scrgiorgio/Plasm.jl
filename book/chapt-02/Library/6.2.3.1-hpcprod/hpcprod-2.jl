function cellprod(cells1, cells2, V, W)
   cells = Vector[]
   for c1 in cells1
      for c2 in cells2
         cell = Vector[]
         for vc in c1
            for wc in c2
               push!(cell, [V[:,vc];W[:,wc]] )
            end
         end
         append!(cells, [cell])
      end
   end
   return cells
end

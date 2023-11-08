function vertprod(V, W)
   vertices = DataStructures.OrderedDict(); k = 0
   for j in 1:size(V,2)
      v = V[:,j]
      for i in 1:size(W,2)
         w = W[:,i]
         id = [v; w]; k += 1
         vertices[id] = k
      end
   end
   return vertices
end

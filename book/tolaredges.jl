using Plasm

function tolaredges(obj::Hpc)::Vector{Vector{Int64}}
   lar=ToLar(obj)
   V  = lar.V
   FV = lar.C[:FV]

   copVV = lar2cop(FV)' * lar2cop(FV)
   predicate(x) = copVV[x...]==2 && x[1]<x[2]
   edges = filter(x -> predicate(x), sort(collect(zip(findnz(copVV)[1:2]...))))

   EV = map(collect, edges)
   return EV
end



# TODO 

## [DONE] (8-12-24) Estendere VIEWCOMPLEX alle celle solide (atomi) quando note e convesse. 
   In particolare con sintassi tipo:
   `VIEWCOMPLEX( lar, show=["CV"], explode=[2,2,2] ) # DA FARE`
   
## [DONE] (8-12-24) `git merge` (branch paoluzzi): 
ho scritto due file:
devi sostituire `SIMPLEX`  (o eliminarlo)  in `fenvs.jl`
e scommentare i test su `SIMPLEX`.  Ovviamente dovresti inserire i due file nel `main` ...

## [DONE] (6-12-24) Ripristinare SPLIT come l'avevo scritta: ritorna `atoms, outer`
L'errore che mi ponevi ieri nel 99,99% dei casi **concreti** non si può presentare.  

## [TODO:: will take time ] (6-12-24) Eliminare vicoli di tipo troppo stretti

In particolare **tutte** le volte che fissi `::Float64` è meglio che metti `::Number`.
Prima ho dovuto capire perché l'esempio `TetGen.jl`  mi dava errore in `MKPOL`.
Motivo:  tu fissavi `Vector{Vector{Int64}}` mentre i dati erano `Vector{Vector{Int32}}` ...





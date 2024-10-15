# TODO

**Eliminare vicoli di tipo troppo stretti**
In particolare **tutte** le volte che fissi `::Float64` è meglio che metti `::Number`.
Prima ho dovuto capire perché l'esempio `TetGen.jl`  mi dava errore in `MKPOL`.
Motivo:  tu fissavi `Vector{Vector{Int64}}` mentre i dati erano `Vector{Vector{Int32}}` ...



# DONE 

- Estendere VIEWCOMPLEX alle celle solide (atomi) quando note e convesse. 
- `git merge` (branch paoluzzi)
- replace `SIMPLEX` 
- Ripristinare SPLIT come l'avevo scritta: ritorna `atoms, outer`









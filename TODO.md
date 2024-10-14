# TODO

### (14-08-24) Va in errore ...

julia> X = GRID([2.4,4.5,-3,4.5,2.4])
julia> Y = GRID([7,5]); Z = GRID([3,3]);
julia> building = X * Y * Z
julia> lar = ARRANGE3D(LAR(building));

julia> atoms,outer = SPLIT(lar)  # io l'avrei fatto solo per CF ... 

ma ho seguito il tuo esempio di SPLIT (non capisco perché sta nel modulo Plasm.jl)


### (11-08-24) Numerare in VIEWCOMPLEX facce, spigoli vertici
visualizzare complessi a celle Lar scrivendo numerazioni di vertice, spigolo, faccia, sia separatamente che tutti insieme.
Per il momento è più urgente modificare con un parametro globale le dimensioni dei testi numerici che orientare la vista verso l'osservatore

### (06-08-24) Eliminare vicoli di tipo troppo stretti ###
In particolare **tutte** le volte che fissi `::Float64` è meglio che metti `::Number`.
Prima ho dovuto capire perché l'esempio `TetGen.jl`  mi dava errore in `MKPOL`.
Motivo:  tu fissavi `Vector{Vector{Int64}}` mentre i dati erano `Vector{Vector{Int32}}` ...



# DONE 

- Estendere VIEWCOMPLEX alle celle solide (atomi) quando note e convesse. 
- `git merge` (branch paoluzzi)
- replace `SIMPLEX` 
- Ripristinare SPLIT come l'avevo scritta: ritorna `atoms, outer`








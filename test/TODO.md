# Piccole cose da fare 
## in ordine di urgenza 

### (8-12-24) Estendere VIEWCOMPLEX alle celle solide (atomi) quando note e convesse. 
   In particolare con sintassi tipo:
   
   `VIEWCOMPLEX( lar, show = ["CV"], explode=[2,2,2] ) # DA FARE`
   
   vedi file **test/TODO**
   
   #### NOTA: 
   
   Non è necessario passare sempre per `ARRAGE3D` quando gli `Hull` siano noti.

### (8-12-24) `git merge` (branch paoluzzi): 

ho scritto due file:

`src/lar/simplexn.jl`
`test/simplexn.jl`

devi sostituire `SIMPLEX`  (o eliminarlo)  in `fenvs.jl`

e scommentare i test su `SIMPLEX`.  Ovviamente dovresti inserire i due file nel `main` ...

### (6-12-24) Ripristinare SPLIT come l'avevo scritta: ritorna `atoms, outer`

L'errore che mi ponevi ieri nel 99,99% dei casi **concreti** non si può presentare.  

Col tuo esempio hai barato.  C'era un terzo atomo (A,B; Outer)

*  fai tutte le diagonali 3D dei bbox degli atomi 
*  calcola il bordo del complesso (prodotto di matrice FC per vettore unitario)
*  fanne il box.
*  non esiste che la diagonale del bbox 3D di un atomo sia più lunga di quella del bordo.  

Nel caso peggiore (CREATO AD HOC) ovvero NON CONCRETO, ce ne possono essere due equali, se proprio fai uno zig-zag in x,y,z.  Ma c'è sempre il bordo algebrico di cui sopra, che avrà anche altri atomi.

In questo caso, praticamente impossibile, si può calcolare l'area dei due candidati ... sono due 2-cicli, quindi si può triangolare e calcolare.

### (6-12-24) Eliminare vicoli di tipo troppo stretti

In particolare **tutte** le volte che fissi `::Float64` è meglio che metti `::Number`.

Prima ho dovuto capire perché l'esempio `TetGen.jl`  mi dava errore in `MKPOL`.
Motivo:  tu fissavi `Vector{Vector{Int64}}` mentre i dati erano `Vector{Vector{Int32}}` ...





# TODO

### (16-08-24) Fare _bene_ SKELETON multidimensionale

<<<<<<< HEAD
Io eviterei di triangolare (nel Viewer) le facce non triangolari. Ad esempio: 
VIEW(sphere = SPHERE(1)([8,16]))

Se guardi l'Hpc, vedi che "hulls="[] è fatto tutto di vettori di 4 punti (che sono sempre complanari)


### (16-08-24) Fare _bene_ SKELETON multidimensionale

Non è difficile. Mi sono accorto chè sbagliato per celle convesse 4D facendo SCHLEGEL3D ...  Quando vuoi, ti spiego cosa fare con Hpc standard (hulls fully-dimensional) e  qhull 
=======
Non è difficile. 
Mi sono accorto chè sbagliato per celle convesse 4D facendo SCHLEGEL3D ...  
Quando vuoi, ti spiego cosa fare con Hpc standard (hulls fully-dimensional) e  qhull
>>>>>>> main

## GLText
modificare con un parametro globale le 
dimensioni dei testi numerici che orientare la vista verso l'osservatore

### (11-08-24) Display numbers
visualizzare complessi a celle Lar scrivendo numerazioni di vertice, spigolo, faccia, sia separatamente che tutti insieme.


# DONE 

- Estendere VIEWCOMPLEX alle celle solide (atomi) quando note e convesse. 
- `git merge` (branch paoluzzi)
- replace `SIMPLEX` 
- Restore SPLIT
- Removed strict type checking
- Fixed:
```julia
        function CROSSPOLYTOPE(D)
                points = Vector{Number}[]
                for i in 1:D
                        point_pos = [0 for x in 1:D];point_pos[i] = +1
                        point_neg = [0 for x in 1:D];point_neg[i] = -1
                        push!(points, point_pos, point_neg)
                end
                cells = [collect(1:D*2)]
                pols = [ [1] ]
                return MKPOL(points, cells, pols)
        end
        CROSSPOLYTOPE(2)
```
- Fixed:
```julia
        X = GRID([2.4,4.5,-3,4.5,2.4])
        Y = GRID([7,5])
        Z = GRID([3,3]);
        building = X * Y * Z
        lar = ARRANGE3D(LAR(building))
        atoms,outer = SPLIT(lar)  
```







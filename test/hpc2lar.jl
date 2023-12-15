using Plasm

# SPHERE: ::Hpc -> ::Lar
obj1 = SPHERE(1)([3,3]).childs[1].childs[1]
points = obj1.points;
hulls  = obj1.hulls;
facets = obj1.facets;
println("points, hulls, facets: $(LEN(points)), $(LEN(hulls)), $(LEN(facets))")

obj = Lar(SPHERE(1.)([3,3]));
V, FV, EV = obj.V, obj.C[:FV], obj.C[:EV];
println("V, FV, EV: $(size(V,2)), $(LEN(FV)), $(LEN(EV))")

# SPHERE: ::Lar -> ::Hpc
obj2 = ToHPC(V,EV);
VIEW(obj2)
obj3 = ToHPC(V,FV);
VIEW(obj3)

# SPHERE: ::Hpc -> ::Lar
obj4 = Lar(obj3);
V, FV, EV = obj4.V, obj4.C[:FV], obj4.C[:EV];
println("V, FV, EV: $(size(V,2)), $(LEN(FV)), $(LEN(EV))")

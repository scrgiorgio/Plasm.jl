
using Plasm

# SPHERE: ::Hpc -> ::Lar
obj = Lar(SPHERE(1.)([32,64]));
V, FV, EV = obj.V, obj.C[:FV], obj.C[:EV];
println("V, FV, EV: $(size(obj.V,2)), $(LEN(obj.C[:FV])), $(LEN(obj.C[:EV]))");

# SPHERE: ::Lar -> ::Hpc
obj2::Hpc = STRUCT(AA(MKPOL)(DISTL([V[:,k] for k=1:size(V,2)], AA(LIST)(FV))));
VIEW(obj2)

# SPHERE: ::Hpc -> ::Lar
obj3 = Lar(obj2);
V, FV, EV = obj3.V, obj3.C[:FV], obj3.C[:EV];
println("V, FV, EV: $(size(obj.V,2)), $(LEN(obj.C[:FV])), $(LEN(obj.C[:EV]))");

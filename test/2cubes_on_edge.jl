using Plasm
using SparseArrays

rV = [0.70710678118655 0.70710678118655 0.0; 2.829225697485796e-17 0.0 0.0; -0.70710678118655 0.70710678118655 0.0; 0.0 1.4142135623731 0.0; 0.4142135623731001 1.0 0.0; 1.0683217716804808e-17 1.0 0.0; 0.70710678118655 0.70710678118655 1.0; 0.0 0.0 1.0; 0.0 1.4142135623731 1.0; 0.4142135623731001 1.0 1.0; -0.70710678118655 0.70710678118655 1.0; 0.0 1.0 1.0; 1.0 0.0 0.0; 1.0 1.0 0.0; 1.0 0.0 1.0; 1.0 1.0 1.0]
rcopEV = sparse([1, 5, 10, 1, 2, 6, 9, 20, 2, 3, 16, 3, 4, 13, 4, 5, 7, 14, 21, 6, 7, 27, 8, 10, 11, 8, 9, 15, 18, 23, 12, 13, 17, 11, 12, 14, 19, 28, 15, 16, 17, 18, 19, 27, 20, 22, 24, 21, 22, 26, 23, 24, 25, 25, 26, 28], [1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16], Int8[-1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1], 28, 16)
rcopFE = sparse([2, 3, 10, 1, 6, 1, 7, 1, 4, 2, 5, 10, 1, 2, 13, 1, 2, 15, 3, 8, 16, 3, 6, 11, 13, 3, 5, 5, 8, 16, 4, 9, 4, 7, 4, 5, 14, 15, 6, 9, 6, 7, 7, 9, 8, 9, 13, 8, 9, 15, 10, 11, 10, 14, 10, 12, 11, 16, 11, 12, 12, 16, 12, 14, 13, 15, 14, 16], [1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 10, 10, 11, 11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28], Int8[-1, 1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, -1, 1, 1, -1], 16, 28)

FE = cop2lar(rcopFE);
#@show FE;
EV = cop2lar(rcopEV);
#@show EV;
FV = cop2lar(lar2cop(FE) * lar2cop(EV) .÷ Int8(2)); # problem: be careful .. !!
#@show FV;
V = permutedims(rV);
#@show V;
twocubes = MKPOL(V,FV,EV)


model = LAR(twocubes)
VIEWCOMPLEX(model)

FV = [union(CAT([CAT([EV[e]]) for e in f])) for f in FE] ## correct !!
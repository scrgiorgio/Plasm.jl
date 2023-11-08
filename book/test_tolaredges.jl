using Plasm
include("./arrangement.jl")
include("./tolaredges.jl")

twocubes = STRUCT([ CUBE(1), T([1,2,3])([.5,.5,.5]), CUBE(1)])

V  = ToLAR(twocubes).childs[1].points
FV = ToLAR(twocubes).childs[1].facets
CV = ToLAR(twocubes).childs[1].hulls
EV = tolaredges(twocubes::Hpc)

Plasm.VIEW( MKPOL(V, CV) ); 
Plasm.VIEW( MKPOL(V, FV) ); 
Plasm.VIEW( MKPOL(V, EV) );  # TODO: 1. change line color; 

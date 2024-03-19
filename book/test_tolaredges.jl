using Plasm
include("./arrangement.jl")
include("./tolaredges.jl")

twocubes = STRUCT([ CUBE(1), T([1,2,3])([.5,.5,.5]), CUBE(1)])

geo=ToGeometry(twocubes)
V  = geo.points
FV = geo.facets
CV = geo.hulls
EV = tolaredges(twocubes::Hpc)

Plasm.VIEW( MKPOL(V, CV) ); 
Plasm.VIEW( MKPOL(V, FV) ); 
Plasm.VIEW( MKPOL(V, EV) );  # TODO: 1. change line color; 

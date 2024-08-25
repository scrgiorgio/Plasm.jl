# ////// Two unit cubes coincident on the z-edge ///////////////////////////////
export SELECTATOMS, DRAWATOMS

# //////////////////////////////////////////////////////////////////////////////
function SELECTATOMS(V,pols) 
   # extract the vertices of each atom
   bboxes = []
   # extract bounding boxes
   for pol in pols
      ev = pol[1] # atom edges (pairs of vertex indices)
      verts = sort(union(CAT(CAT(ev)))) # atom vertices (vertex indices)
      v = V[verts,:] # atom points 
      xbox = bbox(v); # maxmin of first coordinate
      push!(bboxes, xbox)
   end
   # compute 3D cuboidal volumes as (pmin,pmax)
   boxes = AA(collect∘AA(vec))(bboxes)
   diags = [LinearAlgebra.norm(v2-v1) for (v1,v2) in boxes]
   value, position = findmax(diags)
   outerspace = filter(x->x==pols[position],pols)[1]
   atoms = filter(x->x!=pols[position],pols)
   return outerspace,atoms
end

# //////////////////////////////////////////////////////////////////////////////
function DRAWATOMS(V,copEV,copFE,copCF, pols;outer=true)
#@show (V,copEV,copFE,copCF, pols);
   # Lar computation, with associated dictionaries
   EV = AA(sort)([findnz(copEV[k,:])[1] for k=1:copEV.m]); # vertices per edge
   FE = AA(sort)([findnz(copFE[k,:])[1] for k=1:copFE.m]); # edges per face
   FV = AA(sort)([union(CAT([EV[e] for e in f])) for f in FE]); # verts x face
   W = [V[:,k] for k=1:size(V,2)];   
   Wdict = OrderedDict{Any,Int}(zip(W,1:LEN(W)));
   dictEV = Dict{Vector{Int},Int}(collect(zip(EV,1:length(EV))))
   dictFV = Dict{Vector{Int},Int}(collect(zip(FV,1:length(FV))))
   # extraction of atoms from arranged space
   outerspace,atoms = SELECTATOMS(permutedims(V),pols)   
#@show V
   # Visualization of all atoms
   for k=1:length(atoms)
      localV = (V)
      localEV,localFV = atoms[k] 
      localEV = union(CAT(localEV))
      localFV = AA(sort)(localFV)     
      VIEWCOMPLEX2( localV, localEV, localFV, [string(Wdict[a])  for a in W],
         [string(dictEV[a]) for a in localEV],
         [string(dictFV[a]) for a in localFV],
         properties=Properties("background_color"=>Point4d(1.0,1.0,1.0,1.0) ))
   end # Visualization of outer space (B-rep)
   if outer==true
      localV=(V)
      localEV,localFV = outerspace 
      localEV = union(CAT(localEV))
      localFV = AA(sort)(localFV)   
      VIEWCOMPLEX2( localV, localEV, localFV, [string(Wdict[a])  for a in W],
         [string(dictEV[a]) for a in localEV],
         [string(dictFV[a]) for a in localFV],
         properties=Properties("background_color"=>Point4d(1.0,1.0,1.0,1.0) ))
   end
end

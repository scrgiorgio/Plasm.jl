using Triangulate

export Lar, 
   LAR, 
   HPC,
   LAR_SIMPLEX, 
   LAR_CUBOIDGRID, 
   SKELETON, 
   TRIANGULATE2D, 
   LAR2TRIANGLES,  
   SELECTATOMS, 
   VIEWCOMPLEX, 
   DRAWATOMS

# //////////////////////////////////////////////////////////////////////////////////////
# Linear Algebraic Representation . Data type for Cellular and Chain Complex.
mutable struct Lar
  d::Int # intrinsic dimension
  m::Int # embedding dimension (rows of V)
  n::Int # number of vertices  (columns of V)
  V::Matrix{Float64} # object geometry
  C::Dict{Symbol, AbstractArray} # object topology (C for cells)

  # inner constructors
  Lar() = new( -1, 0, 0, Matrix{Float64}(undef,0,0), Dict{Symbol, AbstractArray}() )
  Lar(m::Int,n::Int) = new( m,m,n, Matrix(undef,m,n), Dict{Symbol,AbstractArray}() )
  Lar(d::Int,m::Int,n::Int) = new( d,m,n, Matrix(undef,m,n), Dict{Symbol,AbstractArray}() ) 
  Lar(V::Matrix) = begin m, n = size(V); new( m,m,n, V, Dict{Symbol,AbstractArray}() ) end
  Lar(V::Matrix,C::Dict) = begin m,n = size(V); new( m,m,n, V, C )  end
  Lar(d::Int,V::Matrix,C::Dict) = begin m,n = size(V); new( d,m,n, V, C )  end
  Lar(d,m,n, V,C) = new( d,m,n, V,C )
end



# //////////////////////////////////////////////////////////////////////////////
# from Hpc -> Lar (trying to keep as many information as possible, still unifying to a unique geometry)
function LAR(obj::Hpc; precision=DEFAULT_PRECISION)::Lar
	geo=ToGeometry(obj, precision=precision)
	n=length(geo.points)    # number of vertices  (columns of V)
	m=length(geo.points[1]) # embedding dimension (rows of V) i.e. number of coordinates
	ret=Lar()
	ret.d=m
	ret.n=n
	ret.m=m 
	ret.V=hcat(geo.points...)
	ret.C[:EV]=geo.edges
	ret.C[:FV]=geo.faces
	ret.C[:CV]=geo.hulls
	return ret
end


# //////////////////////////////////////////////////////////////////////////////
#NOTE: better use MKPOLS to specify what Hpc you want to build 
function HPC(lar::Lar)::Hpc 

	if :FV in keys(lar.C) && length(lar.C[:FV])
		return MKPOLS(lar.V, lar.C[:FV])

	elseif :EV in keys(lar.C) && length(lar.C[:EV])
		return MKPOLS(lar.V, lar.C[:EV])

	else
		error("Empty Lar")
	end

end


# //////////////////////////////////////////////////////////////////////////////
""" return ordered vertices  and edges of the 1-cycle f """
function __find_cycle( EV, FE, f::Int )
   vpairs = [EV[e] for e in FE[f]]
   ordered = []
   (A, B), todo = vpairs[1], vpairs[2:end]
   push!(ordered, A)
   while length(todo) > 0
       found = false
       for (I, (a, b)) in enumerate(todo)
           if a == B || b == B
               push!(ordered, B)
               B = (b == B) ? a : b
               found = true
               deleteat!(todo, I)
               break
           end
       end
       @assert found
   end
   push!(ordered, ordered[1])
   edges = [[a,b] for (a,b) in zip(ordered[1:end-1],ordered[2:end])]
   return Array{Int}(ordered[1:end-1]), edges
end

# //////////////////////////////////////////////////////////////////////////////
""" input old LAR consistent data; output triangulated_faces """
function LAR2TRIANGLES(V, EV, FV, FE) 
   V = size(V,1)==3 ? permutedims(V) : V
   triangulated_faces = Vector{Any}(undef, length(FE))

   for f in 1:length(FE)
      edges_idxs = FE[f]
      edge_num = length(edges_idxs)
      fv, edges = __find_cycle(EV, FE, f) 
      # look for independent vector triple
      points = V[fv,:]
      vmap = Dict(zip(fv,1:length(fv))) # vertex map
      mapv = Dict(zip(1:length(fv),fv)) # inverse vertex map
      edges = [[vmap[A],vmap[B]] for (A,B) in edges]
      v1 = LinearAlgebra.normalize(points[2,:] - points[1,:])
      v2 = [0, 0, 0];   v3 = [0, 0, 0]
      err = 1e-8; i = 3 
      while -err < LinearAlgebra.norm(v3) < err
         v2 = LinearAlgebra.normalize(points[i,:] - points[1,:])
         v3 = LinearAlgebra.cross(v1, v2)
         i = i % size(points,1) + 1
      end
      # independent vector triple in face f 
      M = [v1 v2 v3] 
      projected = (points * M)[:,1:2]
      trias = constrained_triangulation2D(permutedims(projected),edges)  # single face f
      triangulated_faces[f] = [[mapv[v] for v in tria] for tria in trias]

   end
   return triangulated_faces
end


# //////////////////////////////////////////////////////////////////////////////
function LAR_CUBOIDGRID(shape::Vector{Int})::Lar
   obj = INSL(POWER)(AA(GRID1)(shape))
   if RN(obj) == 2
      geo=ToGeometry(obj)
      V  = geo.points
      FV = geo.hulls
      EV = geo.edges
      return Lar(hcat(V...), Dict(:FV=>FV, :EV=>EV))
   else 
      return LAR(obj)
   end
end


# //////////////////////////////////////////////////////////////////////////////
"""
    SKELETON(k::Int)(pol::Hpc)::Hpc
Extract the k-skeleton form `Hpc` value `pol`.
```
"""
function SKELETON(ord::Int)
   function SKELETON0(pol::Hpc)
      larpol = LAR(pol)
      if ord==1
         return MKPOLS(larpol.V, larpol.C[:EV])
      elseif ord==2
         return MKPOLS(larpol.V, larpol.C[:FV])
      elseif ord==3
         return MKPOLS(larpol.V, larpol.C[:CV])
      else error("SKELETON($(ord)) not yet implemented")
      end 
   end
   return SKELETON0
end

# //////////////////////////////////////////////////////////////////////////////
function __simplexfacets(simplices)
  @assert hcat(simplices...) isa Matrix
  out = Array{Int64,1}[]
 for it in simplices
   for v in it
     facet = setdiff(it,v)
     push!(out, facet)
   end
 end
 # remove duplicate faces
 return sort(union(out))
end

# ////////////////////////////////////////////////////////////////
function LAR_SIMPLEX(d; complex=false)
   V = [zeros(d,1) I]
   CV = [collect(1:d+1)]
   C = Dict(Symbol("C$(d)V") => CV)
   if complex == false return Lar(V,C)
   else
      cells = CV
      for k=d:-1:1
      	 cells = __simplexfacets(cells)
	 	 key = Symbol("C$(k-1)V")
         push!(C, (key => cells) )
      end
      return Lar(V,C)
   end
end  

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

# ///////////////////////////////////////////////////
function VIEWCOMPLEX(V , EV,  FV; 
   V_text::Vector{String}=nothing, 
   EV_text::Vector{String}=nothing,
   FV_text::Vector{String}=nothing, 
   properties::Properties=Properties())

  properties["background_color"] = get(properties,"background_color", DEFAULT_LAR_BACKGROUND_COLOR)
  properties["line_color"]       = get(properties,"line_color"      , DEFAULT_LINE_COLOR)
  properties["line_width"]       = get(properties,"line_width"      , DEFAULT_LINE_WIDTH)

  properties["text_v_color" ]    = get(properties,"text_v_color"    , DEFAULT_TEXT_V_COLOR)
  properties["text_ev_color"]    = get(properties,"text_ev_color"   , DEFAULT_TEXT_EV_COLOR)
  properties["text_fv_color"]    = get(properties,"text_fv_color"   , DEFAULT_TEXT_FV_COLOR)

  # font sizes
  properties["v_fontsize"]       = get(properties,"v_fontsize"      , DEFAULT_V_FONTSIZE)
  properties["ev_fontsize"]      = get(properties,"ev_fontsize"     , DEFAULT_EV_FONTSIZE)
  properties["fv_fontsize"]      = get(properties,"fv_fontsize"     , DEFAULT_FV_FONTSIZE)

  if isnothing( V_text)  V_text=[string(I) for I in eachindex( V)] end
  if isnothing(EV_text) EV_text=[string(I) for I in eachindex(EV)] end
  if isnothing(FV_text) FV_text=[string(I) for I in eachindex(FV)] end

  used_vertices=[]
  for v in EV append!(used_vertices,v) end
  for v in FV append!(used_vertices,v) end

  obj = MKPOLS(V,EV)
  batches=Vector{GLBatch}()
  append!(batches,GetBatchesForHpc(obj))

  W = [V[:,k] for k=1:size(V,2)]

  # show vertices
  if properties["text_v_color"][4]>0.0 && properties["v_fontsize"]>0
     for I in eachindex(W)
        append!(batches, GLText( 
           (I in used_vertices ? V_text[I] : ""),
           center=ComputeCentroid([W[it] for it in [I]]), 
           color=properties["text_v_color" ],
           fontsize=properties["v_fontsize"]) )
     end
  end

  # show edges
  if !isnothing(EV) && properties["text_ev_color"][4]>0.0 && properties["ev_fontsize"]>0
     for I in eachindex(EV)
        append!(batches, GLText(EV_text[I],
        center=ComputeCentroid([W[it] for it in EV[I]]), 
        color=properties["text_ev_color"],
        fontsize=properties["ev_fontsize"]) )
     end
  end

  # show faces
  if !isnothing(FV) && properties["text_fv_color"][4]>0.0 && properties["fv_fontsize"]>0
     for I in eachindex(FV)
        append!(batches,
        GLText(FV_text[I],
        center=ComputeCentroid([W[it] for it in FV[I]]), 
        color=properties["text_fv_color"],
        fontsize=properties["fv_fontsize"]))
     end
  end

  View(batches, properties)

end

# //////////////////////////////////////////////////////////////////////////////
function VIEWCOMPLEX(mesh::Lar; properties::Properties=Properties())
  return VIEWCOMPLEX(
     mesh.V,
     :EV in keys(mesh.C) ? mesh.C[:EV] : nothing,
     :FV in keys(mesh.C) ? mesh.C[:FV] : nothing,
     properties=properties
  )
end

# //////////////////////////////////////////////////////////////////////////////
function DRAWATOMS(V,copEV,copFE,copCF, pols;outer=true) # V by cols
    # Lar computation, with associated dictionaries
    EV = AA(sort)([findnz(copEV[k,:])[1] for k=1:copEV.m]); # vertices per edge
    FE = AA(sort)([findnz(copFE[k,:])[1] for k=1:copFE.m]); # edges per face
    FV = AA(sort)([union(CAT([EV[e] for e in f])) for f in FE]); # verts x face
    W = [V[k,:] for k=1:size(V,1)];   
    Wdict = OrderedDict{Any,Int}(zip(W,1:LEN(W)));
    dictEV = Dict{Vector{Int},Int}(collect(zip(EV,1:length(EV))))
    dictFV = Dict{Vector{Int},Int}(collect(zip(FV,1:length(FV))))
    # extraction of atoms from arranged space
    outerspace,atoms = SELECTATOMS(permutedims(V),pols)   
    # Visualization of all atoms
    for k=1:length(atoms)
       localV = (V)
       localEV,localFV = atoms[k] 
       localEV = union(CAT(localEV))
       localFV = AA(sort)(localFV)     
       VIEWCOMPLEX(localV, localEV, localFV,
          V_text=[string(Wdict[a])  for a in W],
          EV_text=[string(dictEV[a]) for a in localEV],
          FV_text=[string(dictFV[a]) for a in localFV])
    end # Visualization of outer space (B-rep)
    if outer==true
       localV=(V)
       localEV,localFV = outerspace 
       localEV = union(CAT(localEV))
       localFV = AA(sort)(localFV)   
       VIEWCOMPLEX( localV, localEV, localFV,
          V_text=[string(Wdict[a])  for a in W],
          EV_text=[string(dictEV[a]) for a in localEV],
          FV_text=[string(dictFV[a]) for a in localFV])
    end
  end
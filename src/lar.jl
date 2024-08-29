using Triangulate

export Lar, Hpc, LAR, MKPOLS, HPC
export simplex, CUBOIDGRID, GRID1, SKELETON, ViewCuboidGrid, VIEWCOMPLEX, VIEWCOMPLEX,SQRT, VIEWCOMPLEX2
export triangulate2d, lar2triangles
export SELECTATOMS, DRAWATOMS

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
function MKPOLS(V::Vector{Vector{Float64}}, hulls::Vector{Vector{Int}})::Hpc  
	out = STRUCT(AA(MKPOL)(DISTL(V, AA(LIST)(hulls)))) 
	return out 
end

function MKPOLS(V::Matrix{Float64}, hulls::Vector{Vector{Int}}) 
	W=[V[:,k] for k=1:size(V,2)]
	return MKPOLS(W, hulls)
end

function MKPOLS(V::Union{Vector{Vector{Float64}}, Matrix{Float64}}, cells::Dict{Symbol, AbstractArray}) 
	v=[]
	for (__symbol, hulls) in cells
		push!(v,MKPOLS(V,hulls))
	end
	return STRUCT(v)
end

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
function lar2triangles(V, EV, FV, FE) 
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
      trias = triangulate2d(permutedims(projected),edges)  # single face f
      triangulated_faces[f] = [[mapv[v] for v in tria] for tria in trias]

   end
   return triangulated_faces
end


# //////////////////////////////////////////////////////////////////////////////
"""
    CUBOIDGRID(shape::Vector{Int})::Lar
Generate a cuboidal ``d``-complex grid.

The cuboidal grid has dimension ``d ≥ 2`` equal to `LEN(shape)`, andnumber of rows, columns, pages, etc., depending on `shape` ``Vector``.
# Examples
```jldoctest
julia> CUBOIDGRID([1,1])
Lar(2, 2, 4, [1.0 0.0 0.0 1.0; 0.0 0.0 1.0 1.0], Dict{Symbol, AbstractArray}(:FV => [[1, 2, 3, 4]], :EV => [[3, 4], [2, 3], [1, 2], [1, 4]]))

julia> CUBOIDGRID([2,1])
Lar(2, 2, 6, [1.0 0.0 … 2.0 2.0; 0.0 0.0 … 0.0 1.0], Dict{Symbol, AbstractArray}(:FV => [[1, 2, 3, 4], [5, 1, 4, 6]], :EV => [[3, 4], [2, 3], [1, 2], [1, 4], [4, 6], [1, 5], [5, 6]]))

julia> CUBOIDGRID([2,1,1])
Lar(3, 3, 12, [1.0 0.0 … 2.0 2.0; 0.0 0.0 … 0.0 1.0; 0.0 0.0 … 1.0 1.0], Dict{Symbol, AbstractArray}(:CV => [[1, 2, 3, 4, 5, 6, 7, 8], [9, 1, 10, 3, 11, 5, 12, 7]], :FV => [[4, 2, 3, 1], [5, 6, 2, 1], [5, 4, 7, 1], ..., [5, 11, 9, 1], [11, 10, 9, 12], [4, 7, 10, 12], [5, 7, 11, 12]], :EV => [[2, 1], [2, 3], [4, 1], ..., [10, 12], [11, 9], [11, 12]]))

julia> VIEWCOMPLEX(CUBOIDGRID([2,1,1]))
```
"""
function CUBOIDGRID(shape::Vector{Int})::Lar
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


# //////
# scrgiorgio: commented: don't need the QHull dependency!
#using QHull
#export CHULL
#function CHULL(points::Matrix)
#   ch = QHull.chull(points)
#   FV = ch.simplices
#   pairs = [map(sort,[[i,j],[i,k],[j,k]]) for (i,j,k) in FV]
#   EV = sort!(union(vcat(pairs...)))
#   ret = Lar(2, Matrix(points'), Dict(:EV => EV, :FV => FV))
#end

# //////////////////////////////////////////////////////////////////////////////

function __simplexfacets(simplices)
  @assert hcat(simplices...) isa Matrix
  out = Array{Int64,1}[]
 for simplex in simplices
   for v in simplex
     facet = setdiff(simplex,v)
     push!(out, facet)
   end
 end
 # remove duplicate faces
 return sort(union(out))
end


# /////////////////////////////////////////////////////////////////////////////
"""
    simplex(d; complex=false)::Lar
Generator of `Lar` simplex object of dimension `d`.

Simplex object of `Lar` type with arbitrary dimension `d`. It is the convex combination of ``d+1`` affinely independent points.

# Arguments 
- `dim::Integer=1`: the dimensions along which to perform the computation.
-  [`complex=false`]: when `true` the whole `boundary` simplicial complex is generated.
# Examples
```julia> simplex(1)
Lar(1, 1, 2, [0.0 1.0], Dict{Symbol, AbstractArray}(:C1V => [[1, 2]]))

julia> simplex(2)
Lar(2, 2, 3, [0.0 1.0 0.0; 0.0 0.0 1.0], Dict{Symbol, AbstractArray}(:C2V => [[1, 2, 3]]))

julia> simplex(3, complex=true)
Lar(3, 3, 4, [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0], Dict{Symbol, AbstractArray}(:C3V => [[1, 2, 3, 4]], :C0V => [[1], [2], [3], [4]], :C2V => [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]], :C1V => [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]))
```
"""
function simplex(d; complex=false)
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
function VIEWCOMPLEX(V , EV,  FV; V_text=nothing, EV_text=nothing,FV_text=nothing, properties::Properties=Properties())

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
  #@show (V,copEV,copFE,copCF, pols);
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
  #@show V
    # Visualization of all atoms
    for k=1:length(atoms)
       localV = (V)
       localEV,localFV = atoms[k] 
       localEV = union(CAT(localEV))
       localFV = AA(sort)(localFV)     
       VIEWCOMPLEX( localV, localEV, localFV,
          V_TEXT=[string(Wdict[a])  for a in W],
          EV_TEXT=[string(dictEV[a]) for a in localEV],
          FV_TEXT=[string(dictFV[a]) for a in localFV])
    end # Visualization of outer space (B-rep)
    if outer==true
       localV=(V)
       localEV,localFV = outerspace 
       localEV = union(CAT(localEV))
       localFV = AA(sort)(localFV)   
       VIEWCOMPLEX( localV, localEV, localFV,
          V_TEXT=[string(Wdict[a])  for a in W],
          EV_TEXT=[string(dictEV[a]) for a in localEV],
          FV_TEXT=[string(dictFV[a]) for a in localFV])
    end
  end
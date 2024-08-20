using LinearAlgebra
using QHull

export IsPolytope, IsSimplex, simplex, simplexfacets, CHULL, CUBOIDGRID, GRID1, SKELETON, ViewCuboidGrid, SPHERE, VIEWCOMPLEX, TORUS, RING,VIEWCOMPLEX,SQRT

import Base.-  
-(f::Function, g::Function) = (x...) -> f(x...) - g(x...)  

import Base.+  
+(f::Function, g::Function) = (x...) -> f(x...) + g(x...)  

import Base./  
/(f::Function, g::Function) = (x...) -> f(x...) / g(x...)  

import Base.*  
*(f::Function, g::Function) = (x...) -> f(x...) * g(x...)  

import Base.*
*(pol1::Hpc, pol2::Hpc) = Power(pol1, pol2)

import Base.^
^(f1::Function, f2::Function) = (x,y) -> f1(x)^f2(y) 

import Base.sqrt
SQRT(f::Function) = x -> f(x)^(1/2) 

# //////////////////////////////////////////////////////////////////////////////
"""
    IsPolytope
Plasm predicate `Expr -> Bool` in pure FL style.

Polytopes are the generalization of three-dimensional polyhedra to any number of dimensions.
"""
IsPolytope = AND ∘ CONS([  # pure FL style
   ISPOL, 
   EQ ∘ CONS([LEN ∘ S2 ∘ UKPOL, K(1)])
])

# //////////////////////////////////////////////////////////////////////////////
"""
    IsSimplex
Plasm predicate `Expr -> Bool` in pure FL style.

generalization of the notion of a triangle or tetrahedron to arbitrary dimensions.
"""
IsSimplex = AND ∘ CONS([  # pure FL style
   IsPolytope, 
   EQ ∘ CONS([LEN ∘ S1 ∘ UKPOL, RN + K(1)]) 
])

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
      	 cells = simplexfacets(cells)
	 	 key = Symbol("C$(k-1)V")
         push!(C, (key => cells) )
      end
      return Lar(V,C)
   end
end  

# //////////////////////////////////////////////////////////////////////////////
"""
    simplexfacets(simplices::Vector{Vector{Int}})::Vector{Vector{Int64}}
Generate the set of all faces of `simplices` vector.

See also [`simplex`](@ref)

# Examples
```jldoctest
julia> S = simplex(3, complex=true)
Lar(3, 3, 4, [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0], Dict{Symbol, AbstractArray}(:C3V => [[1, 2, 3, 4]], :C0V => [[1], [2], [3], [4]], :C2V => [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]], :C1V => [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]))

julia> simplexfacets([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]])
6-element Vector{Vector{Int64}}:
 [1, 2]
 [1, 3]
 [1, 4]
 [2, 3]
 [2, 4]
 [3, 4]
```
"""
function simplexfacets(simplices)
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

# //////////////////////////////////////////////////////////////////////////////
"""
    CHULL(points::Matrix):Lar
Generate the convex hull of a matrix of points.

Possible more details for function role explanation.
# Examples
```jldoctest
julia> points = rand(25, 3);

julia> obj = CHULL(points);

julia> VIEW(MKPOLS(obj.V, obj.C[:EV]))
```
"""
function CHULL(points::Matrix)
   ch = QHull.chull(points)
   FV = ch.simplices
   pairs = [map(sort,[[i,j],[i,k],[j,k]]) for (i,j,k) in FV]
   EV = sort!(union(vcat(pairs...)))
   ret = Lar(2, Matrix(points'), Dict(:EV => EV, :FV => FV))
end

# //////////////////////////////////////////////////////////////////////////////
"""
    GRID1(n::Int)::Hpc
Generate a 1D object of `Hpc` type with `n` unit segments.
```
"""
GRID1(n) = QUOTE(DIESIS(n)(1.0))

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
    VIEWCOMPLEX(mesh::Lar; background=Point4d(1,1,1,1))
Visualize any `Lar` object contained in `mesh` argument.

Remark: Any `Hpc` object may by visualized with numbered cells of its 2D or boundary mesh. 
# Examples
```jldoctest
	controlpoints = [[[ 0,0,0],[0 ,3  ,4],[0,6,3],[0,10,0]],[[ 3,0,2], [2 ,2.5,5],[3,6,5],[4,8,2]], [[ 6,0,2],[8 ,3 , 5],[7,6,4.5],[6,10,2.5]], [[10,0,0],[11,3,4],[11,6,3],[10,9,0]]]
	domain = Power(INTERVALS(1.0)(10),INTERVALS(1.0)(10))
	mapping = BEZIERSURFACE(controlpoints)
	VIEWCOMPLEX(LAR(MAP(mapping)(domain)))
```
"""


# ///////////////////////////////////////////////////
function VIEWCOMPLEX(
   lar_vertices::Matrix{Float64}, 
   edges::Vector{Vector{Int}},  
   faces::Vector{Vector{Int}};
   properties::Properties=Properties())

   batches=Vector{GLBatch}()

   # draw vertices 
   vertices= [lar_vertices[:,k] for k=1:size(lar_vertices,2)]
   begin
      color = get(properties,"text_v_color", Point4d(0.0,0.0,0.0, 1.0))
      fontsize=get(properties,"v_fontsize", DEFAULT_V_FONTSIZE)
      if color[4]>0.0 && fontsize>0

         used_vertices=[]
         for v in edges append!(used_vertices,v) end
         for v in faces append!(used_vertices,v) end

         text=get(properties,"v_text", [string(I) for I in eachindex(vertices)])
         for I in eachindex(vertices)
            # do not show unused vertices
            if I in used_vertices 
               append!(batches, GLText(text[I], center=vertices[I],  color=color, fontsize=fontsize)) 
            end
         end
      end
   end

   # draw edges
   if !isnothing(edges)
      batch = GLBatch(LINES)
      batch.line_width  = get(properties,"line_width", 2)
      batch.line_color  = get(properties,"line_color", DEFAULT_LINE_COLOR)
      for edge in edges
            p0 = ToVector3(vertices[edge[1]])
            p1 = ToVector3(vertices[edge[2]])
            push!(batch.vertices.vector, p0...)
            push!(batch.vertices.vector, p1...)
      end
      push!(batches,batch)

      text=get(properties,"ev_text", [string(I) for I in eachindex(edges)])
      color = get(properties,"text_ev_color"   , Point4d(0.0,0.0,1.0, 1.0))
      fontsize = get(properties,"ev_fontsize", DEFAULT_EV_FONTSIZE)
      if color[4]>0.0 && fontsize > 0
         for I in eachindex(edges)
            centroid=ComputeCentroid([vertices[it] for it in edges[I]])
            append!(batches, GLText(text[I], center=centroid, color=color, fontsize=fontsize) )
         end
      end
   end

   # draw faces. NOTE that faces could be non-convex (!)
   if !isnothing(faces)
      batch = GLBatch(TRIANGLES)
      batch.face_color  = get(properties,"face_color", DEFAULT_FACE_COLOR)
      for face in faces
         @show face
         for S in 1:(length(face)-2)
            p0 = ToVector3(vertices[face[S+0]])
            p1 = ToVector3(vertices[face[S+1]])
            p2 = ToVector3(vertices[face[S+2]])
            n = ComputeTriangleNormal(p0, p1, p2) # in case there are degenerate triangles...
            push!(batch.vertices.vector, p0...);push!(batch.normals.vector, n...)
            push!(batch.vertices.vector, p1...);push!(batch.normals.vector, n...)
            push!(batch.vertices.vector, p2...);push!(batch.normals.vector, n...)
         end
      end
      push!(batches,batch)

      text=get(properties,"fv_text", [string(I) for I in eachindex(faces)])
      color = get(properties,"text_fv_color", Point4d(0.0,0.2,0.6, 1.0))
      fontsize = get(properties,"fv_fontsize", DEFAULT_FV_FONTSIZE)
      if color[4]>0.0
         for I in eachindex(faces)
            centroid=ComputeCentroid([vertices[it] for it in faces[I]])
            append!(batches,GLText(text[I], center=centroid, color=color, fontsize=fontsize))
         end

      end
   end

   View(batches,properties)
end

# ////////////////////////////////////////////////////////////////////////////////////
function VIEWCOMPLEX(mesh::Lar; properties::Properties=Properties())
   return VIEWCOMPLEX(
      mesh.V,
      (:EV in keys(mesh.C)) ? mesh.C[:EV] : nothing,
      (:FV in keys(mesh.C)) ? mesh.C[:FV] : nothing,
      properties=properties
   )
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
"""
    SPHERE(radius=1.0::Number)(subds=[16,32]::Vector{Int})
Generate a polyhedral approximation of a spherical surface in 3D.
Maximum correct refinement is LAR(SPHERE(2)([73,40]))

# Examples
```jldoctest
julia> VIEW(SPHERE()())

julia> VIEWCOMPLEX(LAR(SPHERE(2)([4,8])))
```
"""
function SPHERE(radius=1.0::Number)
	function SPHERE0(subds=[16,32]::Vector{Int})
		N, M = subds
		domain = T(1,2)(-pi/2, -pi)(Power(INTERVALS(pi)(N), INTERVALS(2*pi)(M)))
		fx = p -> radius * (-cos(p[1])) * sin(p[2])
		fy = p -> radius * cos(p[1]) * cos(p[2])
		fz = p -> radius * sin(p[1])
		return MAP([fx, fy, fz])(domain)
	end
	return SPHERE0
end

# //////////////////////////////////////////////////////////////////////////////
"""
    TORUS(radii::Vector=[1.0,2])(subds::Vector{Int}=[16,32]):Hpc
Generate polyhedral approximations of a torus surface in 3D.

# Examples
```jldoctest
julia> VIEW(TORUS()())

julia> VIEWCOMPLEX(LAR(TORUS([1,2.])([4,8])))
```
"""
function TORUS(radii=[1.0,2]::Vector)
	r1, r2 = radii
	function TORUS0(subds=[16,32]::Vector{Int})
		N, M = subds
		a = 0.5*(r2-r1)
		c = 0.5*(r1+r2)
		domain = Power(INTERVALS(2*pi)(N), INTERVALS(2*pi)(M))
		fx = p -> (c+a*cos(p[2])) * cos(p[1])
		fy = p -> (c+a*cos(p[2])) * sin(p[1])
		fz = p -> a*sin(p[2])
		return MAP([fx, fy, fz])(domain)
	end
	return TORUS0
end


const Cell = Vector{Int}
export Cell

const Cells = Vector{Cell}
export Cells

# //////////////////////////////////////////////////////////////////////////////
# Linear Algebraic Representation . Data type for Cellular and Chain Complex.
mutable struct Lar

	# object geometry
	# always stored by column i.e. a point is a column 
	#     x1 x2 x3 ...
	#     y1 y2 y3 ...
	V::Points

	# object topology (C for cells)
	C::Dict{Symbol,AbstractArray}

	# for rendering text
	text::Dict{Symbol, Dict{Int,String}}

	# constructor
	Lar(V::Matrix{Float64}=Matrix{Float64}(undef, 0, 0), C::Dict=Dict{Symbol,AbstractArray}()) = begin
		new(V, C, Dict{Int,String}())
	end

end
export Lar

"""
Important!
usually we represent a LAR complex using FV and EV
but we cannot reconstruct FE combinatorially because it would fail on non-convex or holed-faces

So we should migrate to a representation using FE and EV

"""

# ////////////////////////////////////////////
function Base.show(io::IO, lar::Lar) 
	println(io, "Lar(")
	println(io,"  [ # total ", size(lar.V,2))
	for (P,point) in enumerate(eachcol(lar.V))
		println(io,"    ",join([string(it) for it in point]," "), P<=(size(lar.V,2)-1) ? ";" : "","# ",P)
	end
	println(io,"  ],")
	println(io,"  Dict(")
	for (K,key) in enumerate(keys(lar.C))
		cells=lar.C[key]
		println(io, "  ", repr(key)," =>", "[ ", "# total ",length(cells))
		for (C,cell) in enumerate(cells)
			println(io, "    ", repr(cell), C<=(length(cells)-1) ? "," : "", "# ", C)
		end
		println(io, "  ", "]", K<=(length(keys(lar.C))-1) ? "," : "")
	end
	println(io, "  ))")
end


# //////////////////////////////////////////////////////////////////////////////
""" remove duplicate int indices """
function SIMPLIFY(cell::Cell)::Cell
	return collect(sort(collect(Set(cell))))
end


""" simplify all lar cells by sorting vertex indices and removing duplicates

it also try to keep :CF for exploded view
"""
function SIMPLIFY(lar::Lar)

	ret=Lar(lar.V,Dict())

	# techically I could remove replicated vertices, not doing it

	for key in keys(lar.C)

		# ____________________________ EV
		if key==:EV
			ret.C[:EV]=[]
			mapping=Dict{Cell,Int}() # from vertex indices to new edge id
			for (E,ev) in enumerate(lar.C[:EV])
				ev=SIMPLIFY(ev)
				if !(ev in keys(mapping))
					push!(ret.C[:EV], ev)
					mapping[ev]=length(ret.C[:EV])
				end
			end

			continue
		end

		# ____________________________  FV
		if key==:FV
			ret.C[:FV]=[]
			mapping=Dict{Cell,Int}() # from vertex indices to new face id
			for (F,fv) in enumerate(lar.C[:FV])
				fv=SIMPLIFY(fv)
				if ! (fv in keys(mapping))
					push!(ret.C[:FV], fv)
					mapping[fv]=length(ret.C[:FV])
				end
			end

			if haskey(lar.C,:CF)
				ret.C[:CF]=[]
				for cf in lar.C[:CF]
					v::Cell=[]
					for face in cf
						fv=SIMPLIFY(lar.C[:FV][face])
						append!(v,mapping[fv])
					end
					push!(ret.C[:CF],SIMPLIFY(v))
				end
			end

			continue

		end

		#already handled in FV case (see code above)
		if key==:CF
			continue 
		end	
		
		# should I support this case?
		@assert(false)

	end

	return ret

end
export SIMPLIFY

# //////////////////////////////////////////////////////////////////////////////
function SELECT(src::Lar, selected_faces::Cell)::Lar

  FV=src.C[:FV]
  EV=src.C[:EV]

	# TODO: I need to keep FE which is needed for rendering

  selected_vertices=Set(CAT([FV[F] for F in selected_faces]))
  selected_edges=[E for (E,(a,b)) in enumerate(EV) if a in selected_vertices && b in selected_vertices]

  ret=Lar(src.V, Dict(
    :FV => [FV[F] for F in selected_faces],
    :EV => [EV[E] for E in selected_edges]
    ))

  # name mapping
  ret.text[:FV]=Dict{Int,String}(I => string(F) for (I,F) in enumerate(selected_faces))
  ret.text[:EV]=Dict{Int,String}(I => string(E) for (I,E) in enumerate(selected_edges))

  return ret
end
export SELECT

# //////////////////////////////////////////////////////////////////////
function TRIANGULATE(V::Points, EV::Cells)::Cells

	triin = Triangulate.TriangulateIO()
	triin.pointlist = V 
	triin.segmentlist = hcat(EV...)

	(triout, __vorout) = Triangulate.triangulate("pQ", triin) 

	ret=Cells()
	for (u, v, w) in eachcol(triout.trianglelist)
		centroid = (V[:, u] + V[ :,v] + V[:,w]) ./ 3
		if classify_point(centroid, BYROW(V), EV)=="p_in"
			push!(ret, [u, v, w])
		end
	end
	return ret
end
export TRIANGULATE


# //////////////////////////////////////////////////////////////////////////////
"""from Hpc -> Lar 
"""
function LAR(obj::Hpc; precision=DEFAULT_PRECISION)::Lar
	geo = ToGeometry(obj, precision=precision)
	ret = Lar()
	ret.V = hcat(geo.points...)
	ret.C[:EV] = geo.edges
	ret.C[:FV] = geo.faces
	ret.C[:CV] = geo.hulls
	# to add fe?
	return ret
end
export LAR

# //////////////////////////////////////////////////////////////////////////////
""" Create an Hpc from Lar 

use MKPOLS to specify what Hpc you want to build (like only edges or only 2d-faces)
note: they input arguments must be convex otherwise it will not work....
"""
function HPC(lar::Lar)::Hpc

	if :FV in keys(lar.C) && length(lar.C[:FV])
		return MKPOLS(lar.V, lar.C[:FV])

	elseif :EV in keys(lar.C) && length(lar.C[:EV])
		return MKPOLS(lar.V, lar.C[:EV])

	else
		error("Empty Lar")
	end

end
export HPC

# //////////////////////////////////////////////////////////////////////////////
function CUBOIDGRID(shape::Vector{Int})::Lar
	obj = INSL(POWER)(AA(GRID1)(shape))
	if RN(obj) == 2
		geo = ToGeometry(obj)
		V,FV,EV = geo.points,geo.hulls,geo.edges
		return Lar(hcat(V...), Dict(:FV => FV, :EV => EV))
	else
		return LAR(obj)
	end
end
export CUBOIDGRID


# ////////////////////////////////////////////////////////////////
""" create LAR_SIMPLEX

see also fenv SIMPLEX function
"""
function LAR_SIMPLEX(d; complex=false)

	function simplexfacets(simplices)
		@assert hcat(simplices...) isa Matrix
		out = Array{Int64,1}[]
		for it in simplices
			for v in it
				facet = setdiff(it, v)
				push!(out, facet)
			end
		end
		# remove duplicate faces
		return sort(union(out))
	end

	V = [zeros(d, 1) I]
	CV = [collect(1:d+1)]
	C = Dict(Symbol("C$(d)V") => CV)
	if complex == false
		return Lar(V, C)
	else
		cells = CV
		for k = d:-1:1
			cells = simplexfacets(cells)
			key = Symbol("C$(k-1)V")
			push!(C, (key => cells))
		end
		return Lar(V, C)
	end
end
export LAR_SIMPLEX




# //////////////////////////////////////////////////////////////////////////////
""" can find multiple cycles, a cycle interrupts when the previous edge does not have a common vertex 

e.g two cicles

ret=[[a,b],[b,c],[c,d], [h,k],[k,h],...]
"""
function find_vcycles(EV::Cells)::Cells

	todo=copy(EV)
	
	ret=Cells()
	push!(ret,todo[1])
	todo=todo[2:end]

	while length(todo)>0

		# try to attach to the last cycle
		found=false
		for (I,(a,b)) in enumerate(todo)
			if a == ret[end][2]
				push!(ret, [a,b])
				deleteat!(todo, I)
				found=true
				break
			elseif b == ret[end][2]
				push!(ret, [b,a])
				deleteat!(todo, I)
				found=true
				break
			end
		end

		# create a new cycle
		if !found
			push!(ret,todo[1])
			todo=todo[2:end]
		end

	end

	@assert(length(ret)==length(EV))

	return ret

end
export find_vcycles


# //////////////////////////////////////////////////////////////////////////////
function BATCHES(lar::Lar; show=["V", "EV", "FV", "V_text", "FV_text", "EV_text"], explode=1.5, user_color=nothing)::Vector{GLBatch}

  V = lar.V
  EV = haskey(lar.C, :EV) ? lar.C[:EV] : nothing
  FV = haskey(lar.C, :FV) ? lar.C[:FV] : nothing

  # want V to be 3 dimensional
  if size(V, 1) == 2
    zeros = Matrix{Int64}([0.0 for I in 1:size(V, 2)][:,:]')
		V=vcat(V,zeros)
  end

  batches = Vector{GLBatch}()

  batch_points    = GLBatch(POINTS)
  batch_lines     = GLBatch(LINES)
  batch_triangles = GLBatch(TRIANGLES)

	push!(batches, batch_points)
	push!(batches, batch_lines)
	push!(batches, batch_triangles)

  batch_points.point_size = 4
  batch_lines.line_width  = 2
	batch_triangles.line_width = 0

	Vtext  = haskey(lar.text, :V ) ? lar.text[:V ] : Dict(I => string(I) for I in eachindex(V))
  
  function render_point(v_index::Int, pos, color, text)

		if !(v_index in [3, 6, 9, 12, 14, 16, 17, 18, 19, 20, 21])
			return
		end

    if "V" in show
      push!(batch_points.vertices.vector, pos...)
      push!(batch_points.colors.vector, color...)
    end
    if "V_text" in show
      append!(batches, GLText(text, center=pos, color=color, fontsize=0.04))
    end
  end

	function do_explode(cell_points)
		ret=copy(cell_points)
		centroid = compute_centroid(ret)
		vt = (centroid .* [explode; explode; explode]) - centroid
		for C in 1:size(ret,2)
			ret[:,C]=ret[:,C] .+ vt
		end
		return ret, compute_centroid(ret)
	end

	
  if !isnothing(EV)
		EVtext = haskey(lar.text, :EV) ? lar.text[:EV] : Dict(I => string(I) for I in eachindex(EV))

    for (E,ev) in enumerate(EV)
      cell_points, centroid = do_explode(V[:, ev])
			color = isnothing(user_color) ? RandomColor(E) :  user_color

			if "EV" in show
				render_point(ev[1], cell_points[:,1], color, Vtext[ev[1]]); push!(batch_lines.vertices.vector, cell_points[:,1]...);push!(batch_lines.colors.vector, color...)
				render_point(ev[2], cell_points[:,2], color, Vtext[ev[2]]); push!(batch_lines.vertices.vector, cell_points[:,2]...);push!(batch_lines.colors.vector, color...)
			end
			if "EV_text" in show
				append!(batches, GLText(EVtext[E], center=centroid, color=LIGHT_GRAY, fontsize=0.04))
			end
    end
  end

	
  if !isnothing(FV)

		FVtext = haskey(lar.text, :FV) ? lar.text[:FV] : Dict(I => string(I) for I in eachindex(FV))

    for (F, fv) in enumerate(FV)
      points3d, centroid = do_explode(V[:, fv])
      color = isnothing(user_color) ? RandomColor(F) : user_color

			# how to compute FE if it's not pre-computed? (see Lar note below the class)
			# technically is wrong to compute FE combinatorially
      if "FV" in show && haskey(lar.C,:FE)

				# scrgiorgio: not sure the triangulate supports holes....to verify
				# note that probably after arrangement I am not getting any hole so it should work...
				vmap=Dict(zip(fv,1:length(fv)))

				# WRONG
				# cell_EV = convert(Cells, [[vmap[a],vmap[b]] for (a, b) in EV if a in fv && b in fv])

				cell_EV=Cells()
				fe=lar.C[:FE][F]

				@show(fe)
				@show([EV[E] for E in fe])

				for E in fe
					a,b=EV[E]
					push!(cell_EV,[vmap[a],vmap[b]])
				end
				vcycles = find_vcycles(cell_EV)

				for (v_index, pos) in zip(fv,eachcol(points3d))
					render_point(v_index, pos, color, Vtext[v_index])
				end

				for (a,b) in vcycles
					push!(batch_lines.vertices.vector, points3d[:,a]...);push!(batch_lines.colors.vector, color...)
					push!(batch_lines.vertices.vector, points3d[:,b]...);push!(batch_lines.colors.vector, color...)
				end

        # faces in lar can be anything so I need to triangulate
				begin
					points2d = project_points3d(points3d; double_check=true)(points3d)
					for (u, v, w) in TRIANGULATE(points2d, vcycles)
						p0 = points3d[:,u]
						p1 = points3d[:,v]
						p2 = points3d[:,w]
						n = ComputeTriangleNormal(p0, p1, p2)
						push!(batch_triangles.vertices.vector, p0...);push!(batch_triangles.normals.vector, n...);push!(batch_triangles.colors.vector, color...)
						push!(batch_triangles.vertices.vector, p1...);push!(batch_triangles.normals.vector, n...);push!(batch_triangles.colors.vector, color...)
						push!(batch_triangles.vertices.vector, p2...);push!(batch_triangles.normals.vector, n...);push!(batch_triangles.colors.vector, color...)
					end
				end

      end

      if "FV_text" in show
        append!(batches, GLText(FVtext[F], center=centroid, color=color, fontsize=0.04))
      end

    end
  end

  return batches

end
export BATCHES

# //////////////////////////////////////////////////////////////////////////////
function BATCHES(lars::Vector{Lar}; show=["V", "EV", "FV", "V_text", "FV_text", "EV_text"], explode=1.5)::Vector{GLBatch}
	batches=Vector{GLBatch}()
	for (I,lar) in enumerate(lars)
		append!(batches,BATCHES(lar, show=show, explode=explode, user_color=RandomColor(I)))
	end
	return batches
end


# //////////////////////////////////////////////////////////////////////////////
function RandomLine(size_min::Float64,size_max::Float64)
	size = size_min+rand()*(size_max-size_min)
	return STRUCT(
		T(1,2)(rand(2)...), 
		S([1,2])([size,size]), 
		R([1,2])(2*pi*rand()),
		Plasm.SQUARE(1)
	)
end
export RandomLine

# //////////////////////////////////////////////////////////////////////////////
function RandomBubble()
  vs = rand()
  vt = rand(2)
  return STRUCT(
    T(1,2)(vt...),
    S([1,2])([0.25*vs, 0.25*vs]), 
    CIRCUMFERENCE(1)(rand(3:32))
  )
end
export RandomBubble

# //////////////////////////////////////////////////////////////////////////////
function RandomCube(size_min::Float64,size_max::Float64)
  size = size_min+rand()*(size_max-size_min)
  return STRUCT(
    T(1,2,3)(rand(3)...), 
    S([1,2,3])([size,size,size]), 
    R([1,2])(2*pi*rand()),
    R([2,3])(2*pi*rand()),
    R([1,3])(2*pi*rand()),
    Plasm.CUBE(1) 
  )
end
export RandomCube

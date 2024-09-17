export VIEWCOMPLEX

const Cell = Vector{Int}
export Cell

const Cells = Vector{Cell}
export Cells

# //////////////////////////////////////////////////////////////////////////////
# Linear Algebraic Representation . Data type for Cellular and Chain Complex.
"""
Important!
usually we represent a LAR complex using :FV and :EV
but we cannot reconstruct FE combinatorially because it would fail on non-convex or holed-faces
So we should migrate to a representation using FE and EV
"""
mutable struct Lar

	# object geometry
	# always stored by column i.e. a point is a column 
	#     x1 x2 x3 ...
	#     y1 y2 y3 ...
	V::Points

	# object topology (C for cells)
	C::Dict{Symbol,AbstractArray}

	# for mapping new cell ids to old cell idx
	mapping::Dict{Symbol, Dict{Int,Int}}

	# constructor
	Lar(V::Matrix{Float64}=Matrix{Float64}(undef, 0, 0), C::Dict=Dict{Symbol,AbstractArray}()) = begin
		new(V, C, Dict{Int,Int}())
	end

end
export Lar



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
function lar_used_vertices(lar::Lar)
	ret=[]
	for (a,b) in lar.C[:EV]
		append!(ret,[a,b])
	end
	return remove_duplicates(ret)
end
export lar_used_vertices


# //////////////////////////////////////////////////////////////////////////////
function lar_bounding_box(lar::Lar; only_used_vertices=false)
  V=only_used_vertices ? lar.V[:, CAT(lar.C[:EV]) ] : lar.V
  return collect([vec(it) for it in bbox_create(V)])
end
export lar_bounding_box

# //////////////////////////////////////////////////////////////////////////////
""" remove duplicate int indices """
function remove_duplicates(cell::AbstractVector)::AbstractVector
	# note I want the resuned value to be sorted too (needed for edge for example)
	return collect(sort(collect(Set(cell))))
end
export remove_duplicates


# //////////////////////////////////////////////////////////////////////////////
"""can be used also to simplify

NOTE: ignoring :CF that is generally used to crate the `sel`
"""
function SELECT(lar::Lar, sel::Cell)::Lar

	sel=remove_duplicates(sel)

	ret=Lar(lar.V, Dict(:EV => Cells(), :FV => Cells(), :FE => Cells() ))
	ret.mapping=Dict(:F => Dict{Int,Int}(), :E => Dict{Int,Int}())

	fmap=Set{Vector{Int}}() # from (a,b,c,d,...) 
	emap=Dict{Cell,Int}() # from (a,b) to new edge index
	for Fold in sel

		fv=remove_duplicates(lar.C[:FV][Fold])
		if fv in fmap  continue end

		# add new face
		push!(fmap,fv)
		push!(ret.C[:FV], [])
		push!(ret.C[:FE], [])
		Fnew=length(ret.C[:FV])
		ret.mapping[:F][Fnew]=Fold

		# add FE and FV
		# I need FE here, because I cannot know FE from FV,EV especially for non convex-faces
		@assert(haskey(lar.C,:FE))
		begin
			
			for Eold in lar.C[:FE][Fold]

				(a,b)=remove_duplicates(lar.C[:EV][Eold])

				# already added
				if haskey(emap,[a,b])
					Enew=emap[ [a,b] ]

				# add edge
				else
					push!(ret.C[:EV], [a,b])
					Enew=length(ret.C[:EV])
					emap[ [a,b] ]=Enew
					ret.mapping[:E][Enew]=Eold	
				end

				# adding anyway then I will simplify
				push!(ret.C[:FE][Fnew], Enew)

				push!(ret.C[:FV][Fnew], a)
				push!(ret.C[:FV][Fnew], b)

			end

			ret.C[:FV][Fnew]=remove_duplicates(ret.C[:FV][Fnew])
			ret.C[:FE][Fnew]=remove_duplicates(ret.C[:FE][Fnew])

			@assert(ret.C[:FV][Fnew]==fv)

		end

	end

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
function compute_VE(lar::Lar)
	VE=[ [] for I in 1:size(lar.V,2)]
	for (E,(a,b)) in enumerate(lar.C[:EV])
		if !(E in VE[a]) push!(VE[a],E) end
		if !(E in VE[b]) push!(VE[b],E) end
	end
	return VE
end
export compute_VE

# //////////////////////////////////////////////////////////////////////////////
function compute_FE(lar::Lar; is_convex=False)

	# only convex cells supported
	@assert(is_convex && !haskey(lar.C,:FE)) 

	VE=compute_VE(lar)

	ret = Cells()
  for fv in lar.C[:FV]
    fe=Cell()
    for vertex_index in fv
      for edge_index in VE[vertex_index]
				# if both vertices are in the face, it is an edge of the complex (but this works only if the lar complex is made of complex cells)
        a,b=lar.C[:EV][edge_index]
        if (a in fv) && (b in fv)
          push!(fe,edge_index)
        end
      end
    end
    push!(ret,remove_duplicates(fe))
  end

	return ret
end
export compute_FE

# //////////////////////////////////////////////////////////////////////////////
function compute_FV(lar::Lar)
	@assert(!haskey(lar.C,:FV)) 
	lar.C[:FV]=Cells()
	for (F,fe) in enumerate(lar.C[:FE])
		v=Cell()
		for E in fe append!(v,lar.C[:EV][E]) end
		v=remove_duplicates(v)
		push!(lar.C[:FV],v)
	end
end

# //////////////////////////////////////////////////////////////////////////////
"""from Hpc -> Lar """
function LAR(obj::Hpc; precision=TO_GEOMETRY_DEFAULT_PRECISION_DIGITS)::Lar
	geo = ToGeometry(obj, precision=precision)
	ret = Lar()
	ret.V = hcat(geo.points...)
	ret.C[:EV] = geo.edges
	ret.C[:FV] = geo.faces
	ret.C[:CV] = geo.hulls
	ret.C[:FE] = compute_FE(ret, is_convex=true)  # I need FE for display
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
function run_lar_viewer(viewer::Viewer)
	run_viewer(viewer, properties=Properties(
		"background_color" => Point4d([0.9,0.9,0.9,1.0]),
		"use_ortho" => true,
		"title" => "LAR"
	))
end
export run_lar_viewer

# /////////////////////////////////////////////////////
function get_explosion_vt(points::Points, explode)
	centroid = compute_centroid(points)
	vt = (centroid .* [explode[1]; explode[2]; explode[3]]) - centroid
	return vt
end

# /////////////////////////////////////////////////////
function do_explode(points::Points, vt)
	ret=copy(points)
	for C in 1:size(ret,2) 
		ret[:,C]=ret[:,C] .+ vt 
	end
	return ret
end

# /////////////////////////////////////////////////////
function render_edge(viewer::Viewer, lar::Lar, E::Int; line_color=BLACK, vt=[0.0,0.0,0.0], show=[])

	ev=lar.C[:EV][E]
	edge_points=do_explode(lar.V[:, ev], vt)

	lines,colors=Vector{Float32}(),Vector{Float32}()
	append!(lines, edge_points[:,1]);append!(colors, line_color)
	append!(lines, edge_points[:,2]);append!(colors, line_color)
	render_lines(viewer, lines, colors=colors, line_width=2)

	if "V_text" in show
		render_text(viewer, string(haskey(lar.mapping,:V) ? lar.mapping[:V][ev[1]] : ev[1]), center=edge_points[:,1], color=DARK_GRAY, fontsize=DEFAULT_LAR_FONT_SIZE)
		render_text(viewer, string(haskey(lar.mapping,:V) ? lar.mapping[:V][ev[2]] : ev[2]), center=edge_points[:,2], color=DARK_GRAY, fontsize=DEFAULT_LAR_FONT_SIZE)
	end

	if "EV_text" in show
		render_text(viewer, string(haskey(lar.mapping, :E) ? lar.mapping[:E][E] : I), center=compute_centroid(edge_points), color=LIGHT_GRAY, fontsize=DEFAULT_LAR_FONT_SIZE)
	end
end


# /////////////////////////////////////////////////////
function render_face(viewer::Viewer, lar::Lar, F::Int; face_color=BLACK, vt=[0.0,0.0,0.0], show=[])

	fv=lar.C[:FV][F]
	face_points=do_explode(lar.V[:, fv], vt)
	
	# generic solution is to use FE to compute constrained triangulation and then triangles
	if haskey(lar.C,:FE)
		vmap=Dict(zip(fv,1:length(fv)))
		cell_EV=Cells()
		fe=lar.C[:FE][F]
		for E in fe
			a,b=lar.C[:EV][E]
			push!(cell_EV,[vmap[a],vmap[b]])
		end
		vcycles = find_vcycles(cell_EV)
		points2d = project_points3d(face_points; double_check=true)(face_points)
		triangulation = TRIANGULATE(points2d, vcycles)

	# is it a simple triangle?
	elseif length(fv)==3
		vcycles=[1,2],[2,3],[3,1]
		triangulation=[[1,2,3]]

	else
		# I need to know FE (which cannot be computed automatically from FV EV considering non-convex faces)
		return 
	end

	# render points
	begin
		for (v_index, pos) in zip(fv,eachcol(face_points))
			if "V_text" in show
				render_text(viewer, string(haskey(lar.mapping,:V) ? lar.mapping[:V][v_index] : v_index), center=pos, color=DARK_GRAY, fontsize=0.04)
			end
		end		
		
	end

	# render lines
	begin
		lines, colors=Vector{Float32}(),Vector{Float32}()
		begin
			for (a,b) in vcycles
				append!(lines, face_points[:,a]);append!(colors, DARK_GRAY)
				append!(lines, face_points[:,b]);append!(colors, DARK_GRAY)
			end			
		end
		render_lines(viewer, lines, colors=colors)
	end

	# render triangles
	begin
		triangles, colors=Vector{Float32}(),Vector{Float32}()
		begin
			for (u, v, w) in triangulation
				p0 = face_points[:,u];append!(triangles, p0);append!(colors, face_color)
				p1 = face_points[:,v];append!(triangles, p1);append!(colors, face_color)
				p2 = face_points[:,w];append!(triangles, p2);append!(colors, face_color)
			end
		end
		render_triangles(viewer, triangles, colors=colors, enable_polygon_offset=true)
	end

	# render text
	begin
		if "FV_text" in show
			centroid=compute_centroid(face_points)
			render_text(batches, string(haskey(lar.mapping, :F) ? lar.mapping[:F][F] : F), center=compute_centroid(centroid), color=face_color, fontsize=0.04)
		end
	end

end


# //////////////////////////////////////////////////////////////////////////////
function render_lar(viewer::Viewer, lar::Lar; show=["V", "EV", "FV"], explode=[1.0,1.0,1.0])

  # want lar to be 3 dimensional 
	begin
		lar = Lar(lar.V, lar.C)
		if size(lar.V, 1) == 2
			zeros = Matrix{Int64}([0.0 for I in 1:size(lar.V, 2)][:,:]')
			lar.V=vcat(lar.V,zeros)
		end
	end

	# ____________________________________________
  if "FV" in show && haskey(lar.C,:FV)

		# explode by cell
		if "atom" in show && haskey(lar.C, :CF)
			for (C,cf) in enumerate(lar.C[:CF])
				v_indices=[]
				for F in cf append!(v_indices,lar.C[:FV][F]) end
				v_indices=remove_duplicates(v_indices)
				vt=get_explosion_vt(lar.V[:, v_indices], explode)
				face_color = RandomColor()
				for F in cf
					render_face(viewer, lar, F, vt=vt, face_color=face_color, show=show)
				end
			end
		# expode by single face
		else
			for (F, fv) in enumerate(lar.C[:FV])
				vt=get_explosion_vt(lar.V[:, fv], explode)
				face_color = RandomColor()
				render_face(viewer, lar, F, vt=vt, face_color=face_color, show=show)
			end
		end

	# show lines
  elseif "EV" in show && haskey(lar.C,:EV)

		# explode by polygon
		if "atom" in show && haskey(lar.C, :FV) && haskey(lar.C, :FE)
			for (F,fv) in enumerate(lar.C[:FV])
				vt=get_explosion_vt(lar.V[:, fv], explode)
				line_color = RandomColor()
				for E in lar.C[:FE][F]
					render_edge(viewer, lar, E, vt=vt, line_color=line_color, show=show)
				end
			end
		# explode by single edge
		else
			for (E,ev) in enumerate(lar.C[:EV])
				vt=get_explosion_vt(lar.V[:, ev], explode)
				line_color = RandomColor()
				render_edge(viewer, lar, E, vt=vt, line_color=line_color, show=show)
			end
		end
	end

end
export render_lar

# //////////////////////////////////////////////////////////////////////////////
function VIEWCOMPLEX(viewer::Viewer, lar::Lar; show=["V", "EV", "FV"], explode=[1.0,1.0,1.0])
	render_lar(viewer, lar, show=show, explode=explode)
	run_lar_viewer(viewer)
end

function VIEWCOMPLEX(lar::Lar; show=["V", "EV", "FV"], explode=[1.0,1.0,1.0])
	VIEWCOMPLEX(Viewer(), lar,show=show, explode=explode)
end

export VIEWCOMPLEX

# //////////////////////////////////////////////////////////////////////////////

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


# /////////////////////////////////////////////////////
# From here the sparse part
# /////////////////////////////////////////////////////


const Chain        = SparseVector{Int8,Int}
const ChainOp      = SparseMatrixCSC{Int8,Int}
const ChainComplex = Vector{ChainOp}

export Chain, ChainOp, ChainComplex

# ///////////////////////////////////////////////////////////
""" converte dense to sparse"""
function lar2cop(cells::Cells)::ChainOp
	I, J, V = Int[], Int[], Int[]
	for (C,cell) in enumerate(cells)
		for K in cell
			push!(I, C)
			push!(J, K)
			push!(V, 1)
		end
	end
	return sparse(I, J, V)
end
export lar2cop


""" converte sparse to dense"""
function cop2lar(cop::ChainOp)::Cells
	return [findnz(cop[k, :])[1] for k = 1:size(cop, 1)]
end
export cop2lar

""" EV dense to sparse """
function cop_coboundary_0(EV::Cells)::ChainOp
	copEV = lar2cop(EV)
	copVE = copEV'
	for (E,ev) in enumerate(EV)
		v1,v2=ev
		copVE[v1, E] = -1 # from +1 -> -1
	end
	copEV=LinearAlgebra.transpose(copVE)
	return convert(ChainOp,copEV)
end
export cop_coboundary_0




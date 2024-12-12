export VIEWCOMPLEX

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
	C::Dict{Symbol,Cells}

	# for mapping new cell ids to old cell idx
	mapping::Dict{Symbol, Dict{Int,Int}}

	# extra properties
	properties::Properties

	# constructor
	Lar(V::Matrix{Float64}=Matrix{Float64}(undef, 0, 0), C::Dict=Dict{Symbol,Cells}(); mapping=Dict{Int,Int}(), properties=Properties()) = begin
		new(V, C, mapping, properties)
	end

end
export Lar

# ////////////////////////////////////////////
function lar_copy(src::Lar)::Lar
	# note: sharing vertices (!)
	return Lar(
		src.V,
		deepcopy(src.C),
		mapping=deepcopy(src.mapping)
	)
end
export lar_copy

# ////////////////////////////////////////////
function lar_vertex_name(lar::Lar, id::Int)::String
	return string(haskey(lar.mapping,:V) ? lar.mapping[:V][id] : id)
end

function lar_edge_name(lar::Lar, id::Int)::String
	return string(haskey(lar.mapping,:E) ? lar.mapping[:E][id] : id)
end 

function lar_face_name(lar::Lar, id::Int)::String
	return string(haskey(lar.mapping,:F) ? lar.mapping[:F][id] : id)
end 

# ////////////////////////////////////////////
# disabled: https://github.com/scrgiorgio/Plasm.jl/issues/17

function show_debug(lar::Lar) 
	println("Lar(")
	println("  [ # total ", size(lar.V,2))
	for (P,point) in enumerate(eachcol(lar.V))
		println("    ",join([string(it) for it in point]," "), P<=(size(lar.V,2)-1) ? ";" : "","# ",P)
	end
	println("  ] ,")
	println("   Dict(")
	for (K,key) in enumerate(keys(lar.C))
		cells=lar.C[key]
		println("  ", repr(key)," =>", "[ ", "# total ",length(cells))
		for (C,cell) in enumerate(cells)
			println("    ", repr(cell), C<=(length(cells)-1) ? "," : "", "# ", C)
		end
		println("  ", "]", K<=(length(keys(lar.C))-1) ? "," : "")
	end
	println("  ),")
	println("   Mapping(")
	println( lar.mapping)
	println("  )")
	println(")")
end
export show_debug

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
	return bbox_create(V)
end
export lar_bounding_box



# ////////////////////////////////////////////////////////
function lar_simplify_internal(dst::Lar, src::Lar, symbol::Symbol; allow_duplicates=false)

  if !haskey(src.C, symbol)
    return
  end

  @assert(!haskey(dst.C, symbol))
  dst.C[symbol] = Cells()

  if symbol == :EV
    @assert(!haskey(src.mapping, :E)) # otherwise I would need mapping of mapping
    dst.mapping[:E] = Dict{Int,Int}()
    Dmap = dst.mapping[:E]
    Smap = nothing

  elseif symbol == :FV
    @assert(!haskey(src.mapping, :F)) # otherwise I would need mapping of mapping
    dst.mapping[:F] = Dict{Int,Int}()
    Dmap = dst.mapping[:F]
    Smap = nothing

  elseif symbol == :FE
    @assert(haskey(dst.mapping,:E))
    Dmap = nothing
    Smap = dst.mapping[:E]

  elseif symbol == :CF
    @assert(haskey(dst.mapping,:F))
    Dmap = nothing
    Smap = dst.mapping[:F]

  elseif symbol == :CV
    Dmap = nothing
    Smap = nothing

  else
    @assert(false)

  end

  added = Dict()
  for (old_id, cell) in enumerate(src.C[symbol])
    cell = Cell([isnothing(Smap) ? it : Smap[it] for it in cell])
		simplified = normalize_cell(cell)

		if allow_duplicates || !haskey(added, simplified)
			push!(dst.C[symbol], simplified)
			added[simplified] = length(dst.C[symbol])
		end

		if !isnothing(Dmap)
			Dmap[old_id] = added[simplified]
		end

  end

end

# ////////////////////////////////////////////////////////
function SIMPLIFY(src::Lar)

  dst = Lar(src.V, Dict{Symbol,Cells}())

  # to I need to support other cases?
  @assert(all([key in [:CV, :EV, :FV, :FE, :CF] for key in keys(src.C)]))

  # important the order (so that E mapping comes before :FE and F mapping comes before :CF)
	lar_simplify_internal(dst, src, :EV, allow_duplicates=false) 
	lar_simplify_internal(dst, src, :FV, allow_duplicates=false) 
	lar_simplify_internal(dst, src, :CV, allow_duplicates=false) 
	lar_simplify_internal(dst, src, :FE, allow_duplicates=false) 
	lar_simplify_internal(dst, src, :CF, allow_duplicates=true)  # i need to keep duplicates for inside/outside

	# do not need to keep the mapping here
	dst.mapping=Dict{Symbol, Dict{Int,Int}}()

  return dst
end
export SIMPLIFY



# /////////////////////////////////////////////////////////////////////
function CHECK(lar::Lar)

  # each vertex should have >=2 edges in 2dim and >=3 edges in 3d
  begin
    pdim=size(lar.V,1)
    count=Dict{Int64, Int64}(P => 0 for P in 1:size(lar.V,2))
    for (a,b) in lar.C[:EV]
      count[a]+=1
      count[b]+=1
    end
    @assert(all([v==0 || v>=pdim for (k,v) in count]))
  end

  # eatch edge should have 2 vertives
  for ev in lar.C[:EV]
    @assert(length(ev)==2)
  end

  # each edge should have at least one face
  begin
    count_faces=Dict{Int,Set{Int}}()
    for (F,fe) in enumerate(lar.C[:FE])
      for E in fe
        if !haskey(count_faces,E) count_faces[E]=Set{Int}() end
        push!(count_faces[E],F)
      end
    end
    @assert(all([length(it)>=1 for it in values(count_faces)]))
  end

  # one face should have at least 3 edges
  begin
    for (F,fe) in enumerate(lar.C[:FE])
      @assert(length(fe)>=3)
      # and should be able to find face cycles
      cell_ev=[lar.C[:EV][E] for E in fe]
			#println("F ",F," fe ",fe, " ",  cell_ev)
      find_vcycles(cell_ev)
    end
  end

end

# ///////////////////////////////////////////////////////////////
function SELECT_FACES(lar::Lar, sel::Cell)::Lar

	sel=normalize_cell(sel)

	ret=Lar(lar.V, Dict(
		:EV => Cells(), 
		:FV => Cells(), 
		:FE => Cells() ))

	ret.mapping=Dict(
		:F => Dict{Int,Int}(), 
		:E => Dict{Int,Int}()
	)

	FV=haskey(lar.C,:FV) ? lar.C[:FV] : compute_FV(lar)

	added_fv=Dict{Cell,Int}() # from (a,b,c,d,...)
	added_ev=Dict{Cell,Int}() # from (a,b) to new edge index

	for Fold in sel

		fv=normalize_cell(FV[Fold])

		# already added
		if haskey(added_fv,fv)  
			continue 
		end

		# add new face
		begin
			added_fv[fv]=Fold
			push!(ret.C[:FV], [])
			push!(ret.C[:FE], [])
			Fnew=length(ret.C[:FV])
			ret.mapping[:F][Fnew]=Fold
		end

		# add FE and FV (NOTE: I need FE here, because I cannot know FE from FV,EV especially for non convex-faces)
		@assert(haskey(lar.C,:FE))
		begin
			
			for Eold in lar.C[:FE][Fold]
				(a,b)=normalize_cell(lar.C[:EV][Eold])
				# already added
				if haskey(added_ev,[a,b])
					Enew=added_ev[ [a,b] ]
				# add edge
				else
					push!(ret.C[:EV], [a,b])
					Enew=length(ret.C[:EV])
					added_ev[ [a,b] ]=Enew
					ret.mapping[:E][Enew]=Eold	
				end
				# adding anyway then I will simplify
				push!(ret.C[:FE][Fnew], Enew)
				push!(ret.C[:FV][Fnew], a)
				push!(ret.C[:FV][Fnew], b)
			end
			ret.C[:FV][Fnew]=normalize_cell(ret.C[:FV][Fnew])
			ret.C[:FE][Fnew]=normalize_cell(ret.C[:FE][Fnew])
			@assert(ret.C[:FV][Fnew]==fv) 
		end
	end

	return ret

end

export SELECT_FACES


# //////////////////////////////////////////////////////////////////////
function TRIANGULATE(V::Points, EV::Cells)::Cells

	triin = Triangulate.TriangulateIO()
	triin.pointlist = V 
	triin.segmentlist = hcat(EV...)

	# p: Triangulates a Planar Straight Line Graph
	# Q for quiet
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
function compute_FE(lar::Lar; is_convex=false)

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
    push!(ret,normalize_cell(fe))
  end

	return ret
end
export compute_FE

# //////////////////////////////////////////////////////////////////////////////
function compute_FV(lar::Lar)::Cells 
	ret=Cells()
	for (F,fe) in enumerate(lar.C[:FE])
		v=Cell()
		for E in fe 
			append!(v,lar.C[:EV][E]) 
		end
		v=normalize_cell(v)
		push!(ret,v)
	end
	return ret
end

# //////////////////////////////////////////////////////////////////////////////
function compute_CF(lar::Lar; is_convex=false)::Cells
	@assert(is_convex)
	@assert(haskey(lar.C,:CV))
	@assert(haskey(lar.C,:FV))

	ret=Cells()
	for (C,cv) in enumerate(lar.C[:CV])
		cf=Cell()
		s_cv=Set(cv)
		for (F,fv) in enumerate(lar.C[:FV])
			s_fv=Set(fv)
			if intersect(s_fv, s_cv)==s_fv # face must be completely in the hull
				push!(cf,F)
			end
		end
		push!(ret, cf)
	end
	return simplify_cells(ret)

end

# //////////////////////////////////////////////////////////////////////////////
function compute_CV(lar::Lar)::Cells
	@assert(haskey(lar.C,:CF))
	@assert(haskey(lar.C,:FV)) # NOTE:I think using FV will work only for convex cells, but not sure
	ret=Cells()
	for (C,cf) in enumerate(lar.C[:CF])
		cv=Cell()
		for (F, fv) in enumerate(lar.C[:FV][C])
			append!(cv, fv)
		end
		push!(ret,cv)
	end
	return simplify_cells(ret)
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
	ret.C[:FE] = compute_FE(ret, is_convex=true) # I need for display
	ret.C[:CF] = compute_CF(ret, is_convex=true) #  I need for display
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
function CUBOIDGRID(shape::Vector)::Lar
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


# //////////////////////////////////////////////////////////////////////////////
""" can find multiple cycles, a cycle interrupts when the previous edge does not have a common vertex 

e.g two cicles

ret=[[a,b],[b,c],[c,d], [h,k],[k,h],...]
"""
function find_vcycles(EV::Cells)::Vector{Cells}

	todo=copy(EV)
	
	ret=Vector{Cells}()
	push!(ret,Cells())
	
	push!(ret[end],todo[1])
	todo=todo[2:end]

	while length(todo)>0

		# try to attach to the last cycle
		found=false
		last=ret[end][end][2]
		for (I,(a,b)) in enumerate(todo)
			if a == last
				push!(ret[end],[a,b])
				deleteat!(todo, I)
				found=true
				break
			elseif b == last
				push!(ret[end],[b,a])
				deleteat!(todo, I)
				found=true
				break
			end
		end

		# create a new cycle
		if !found
			push!(ret,Cells())
			push!(ret[end],todo[1])
			todo=todo[2:end]
		end

	end

	@assert(length(vcat(ret...))==length(EV))
	@assert(all([length(it)>=3 for it in ret]))

	return ret

end
export find_vcycles

# //////////////////////////////////////////////////////////////////////////////
function run_lar_viewer(viewer::Viewer; title="LAR", use_thread=false, properties=nothing)

	if isnothing(properties)
		properties=Properties()
	end

	properties["show_axis"]           = get(properties,"show_axis", true)
	properties["title"]               = get(properties,"title",title)
	properties["background_color"]    = get(properties,"background_color",Point4d([0.9,0.9,0.9,1.0]))

	run_viewer(viewer, 
		properties=properties,
		use_thread=use_thread
	)
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
function render_edge(viewer::Viewer, lar::Lar, E::Int; line_color=BLACK, vt=[0.0,0.0,0.0], show=[], properties=Properties())

	ev=lar.C[:EV][E]
	edge_points=do_explode(lar.V[:, ev], vt)

	lines,colors=Vector{Float32}(),Vector{Float32}()
	append!(lines, edge_points[:,1]);append!(colors, line_color)
	append!(lines, edge_points[:,2]);append!(colors, line_color)
	render_lines(viewer, lines, colors=colors, line_width=DEFAULT_LINE_WIDTH)

	if "Vtext" in show
		render_text(viewer, lar_vertex_name(lar, ev[1]), center=edge_points[:,1], color=LAR_VERTEX_COLOR, fontsize=get(properties,"font_size", LAR_VERTEX_FONT_SIZE))
		render_text(viewer, lar_vertex_name(lar, ev[2]), center=edge_points[:,2], color=LAR_VERTEX_COLOR, fontsize=get(properties,"font_size", LAR_VERTEX_FONT_SIZE))
	end

	if "Etext" in show
		render_text(viewer, lar_edge_name(lar, E), center=compute_centroid(edge_points), color=LAR_EDGE_COLOR, fontsize=get(properties,"font_size", LAR_EDGE_FONT_SIZE))
	end
end


# /////////////////////////////////////////////////////
function render_face(viewer::Viewer, lar::Lar, F::Int; face_color=BLACK, vt=[0.0,0.0,0.0], show=[], properties=Properties())

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
		vcycles = vcat(find_vcycles(cell_EV)...)
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
			if "Vtext" in show
				perturbation= [0.0,0.0,0.0] 
				# perturbation=rand(3)*0.01
				render_text(viewer, lar_vertex_name(lar, v_index), center=pos + perturbation , color=LAR_VERTEX_COLOR, fontsize=get(properties,"font_size", LAR_VERTEX_FONT_SIZE))
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
	if face_color[4]>0.0
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
		if "Ftext" in show
			render_text(viewer, lar_face_name(lar, F), center=compute_centroid(face_points), color=LAR_FACE_COLOR, fontsize=get(properties,"font_size", LAR_FACE_FONT_SIZE))
		end
	end

end


# //////////////////////////////////////////////////////////////////////////////
function render_lar(viewer::Viewer, lar::Lar; show=["V", "EV", "FV"], explode=[1.0,1.0,1.0], face_color=nothing, properties=Properties())

  # want lar to be 3 dimensional 
	begin
		lar = lar_copy(lar)
		if size(lar.V, 1) == 2
			zeros = Matrix{Int64}([0.0 for I in 1:size(lar.V, 2)][:,:]')
			lar.V=vcat(lar.V,zeros)
		end
	end

	EV=haskey(lar.C, :EV) && length(lar.C[:EV])>0 ? lar.C[:EV] : nothing
	FE=haskey(lar.C, :FE) && length(lar.C[:FE])>0 ? lar.C[:FE] : nothing
	FV=haskey(lar.C, :FV) && length(lar.C[:FV])>0 ? lar.C[:FV] : nothing
	CF=haskey(lar.C, :CF) && length(lar.C[:CF])>0 ? lar.C[:CF] : nothing

	# show atoms exploding by cell centroid
	if ("CV" in show || "CF" in show) && !isnothing(CF) && !isnothing(FV) 

		for (C,cf) in enumerate(CF)
			v_indices=[]
			for F in cf 
				append!(v_indices,FV[F]) 
			end
			v_indices=remove_duplicates(v_indices)
			vt=get_explosion_vt(lar.V[:, v_indices], explode)
			atom_face_color = isnothing(face_color) ? RandomColor() : face_color
			for F in cf
				render_face(viewer, lar, F, vt=vt, face_color=atom_face_color, show=show, properties=properties)
			end
		end

	end

	# show faces expoding by face centroid
	if "FV" in show && !isnothing(FV)

		for (F, fv) in enumerate(FV)
			vt=get_explosion_vt(lar.V[:, fv], explode)
			# sometimes getting StackOverflowError
			try
				render_face(viewer, lar, F, vt=vt, face_color=isnothing(face_color) ? RandomColor() : face_color, show=show, properties=properties)
			catch
				println("ERROR in render_face, what is going on???")
			end
		end
	end

	# show edges
	if "EV" in show
		
		# exploding by face centroid
		if !isnothing(FV) !isnothing(FE) 
			for (F,fv) in enumerate(FV)
				vt=get_explosion_vt(lar.V[:, fv], explode)
				# line_color = RandomColor()
				line_color=BLACK
				for E in FE[F]
					render_edge(viewer, lar, E, vt=vt, line_color=line_color, show=show, properties=properties)
				end
			end
			
		# exploded by edge centroid
		elseif !isnothing(EV)
			for (E,ev) in enumerate(EV)
				vt=get_explosion_vt(lar.V[:, ev], explode)
				# line_color = RandomColor()
				line_color=BLACK
				render_edge(viewer, lar, E, vt=vt, line_color=line_color, show=show, properties=properties)
			end
		end

		# show Ftext if needed
		if "Ftext" in show && !isnothing(FV)
			for (F, fv) in enumerate(FV)
				vt=get_explosion_vt(lar.V[:, fv], explode)
				# sometimes getting StackOverflowError
				try
					render_face(viewer, lar, F, vt=vt, face_color=TRANSPARENT, show=show, properties=properties)
				catch
					println("ERROR in render_face, what is going on???")
				end
			end

		end

	end

end
export render_lar

# //////////////////////////////////////////////////////////////////////////////
function VIEWCOMPLEX(viewer::Viewer, lar::Lar; show=["V", "EV", "FV"], explode=[1.0,1.0,1.0], face_color=nothing, title="LAR", use_thread=false, properties=Properties(), numbering=false)

	if numbering
		show=[show ;  ["Vtext", "Etext", "Ftext"] ]
	end

	render_lar(viewer, lar, show=show, explode=explode, face_color=face_color, properties=properties)
	run_lar_viewer(viewer, title=title, use_thread=use_thread, properties=properties)
end

function VIEWCOMPLEX(lar::Lar; show=["V", "EV", "FV"], explode=[1.0,1.0,1.0], face_color=nothing, title="LAR", use_thread=false, properties=Properties(), numbering=false)
	VIEWCOMPLEX(Viewer(), lar, show=show, explode=explode, face_color=face_color, title=title, use_thread=use_thread, properties=properties, numbering=numbering)
end

export VIEWCOMPLEX


# /////////////////////////////////////////////////////////////////////
function VIEWEDGES(V::Points, EV::Cells; explode=[1.0,1.0,1.0], title="LAR")
  lar=Lar(V,Dict{Symbol,Cells}(:EV=>EV))
  # print_matrix("sorted_points", hcat(collect(sort([it for it in eachcol(V)]))))
  VIEWCOMPLEX(lar, explode=explode, show=["V","EV","Vtext"], title=title)
end

# /////////////////////////////////////////////////////////////////////
function VIEWEDGES(V::Points, segmentlist::Matrix; explode=[1.0,1.0,1.0], title="LAR")
  VIEWEDGES(V,[Cell(it) for it in eachcol(segmentlist)], explode=explode, title=title)
end

# /////////////////////////////////////////////////////////////////////
function VIEWTRIANGLES(V::Points, triangles::Cells; explode=[1.0,1.0,1.0], title="LAR")
  lar=Lar(V, Dict{Symbol,Cells}(:EV => Cells(),:FE => Cells()))
  for (u,v,w) in triangles
    E=length(lar.C[:EV])
    append!(lar.C[:EV], [[u,v],[v,w],[w,u]])
    push!(lar.C[:FE], [E+1,E+2,E+3])
  end
  lar.C[:FV]=compute_FV(lar)
  VIEWCOMPLEX(lar, explode=explode, show=["V", "EV", "FV", "Vtext"], title=title)
end

function VIEWTRIANGLES(V::Points, triangles::Matrix; explode=[1.0,1.0,1.0], title="LAR")
  return VIEWTRIANGLES(V,[Cell(it) for it in eachcol(triangles)], explode=explode, title=title)
end

# //////////////////////////////////////////////////////////////
function VIEWEDGE(lar::Lar, ev::Cell; enlarge=nothing)

  a,b=ev
  @show(lar)
  for (I,fv) in enumerate(lar.C[:FV])
    println("fv ", fv, " # ",I)
  end
  
  selected_vertices=[a,b]

  # take any vertex near by
  if !isnothing(enlarge)
    for (P,p) in enumerate(eachcol(lar.V))
      if LinearAlgebra.norm(p - lar.V[:,a])<enlarge || LinearAlgebra.norm(p - lar.V[:,b])<enlarge push!(selected_vertices,P) end
    end
  end

  # show all faces touching the vertices
  selected_faces=Cell()
  for (F,fv) in enumerate(lar.C[:FV])
    inters=intersect(Set(selected_vertices),Set(fv))
    if length(inters)>0
      push!(selected_faces,F)
    end
  end

  selected_faces=normalize_cell(selected_faces)
  VIEWCOMPLEX(SELECT_FACES(lar, selected_faces), show=["FV", "Vtext"], explode=[1.0,1.0,1.0], face_color=TRANSPARENT, title="debug edge")

end


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




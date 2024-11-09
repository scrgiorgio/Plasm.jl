




#//////////////////////////////////////////////////////////////////////////////
# REMARK: will not work when vs[:,1:3] are aligned !!!!  TODO: fix 
function fragment_submanifold_mapping(vs)
	""" Compute the map from vs (at least three) to z=0 """
	u1 = vs[2, :] - vs[1, :]
	u2 = vs[3, :] - vs[1, :]
	u3 = LinearAlgebra.cross(u1, u2)
	T = Matrix{Float64}(LinearAlgebra.I, 4, 4)
	T[4, 1:3] = -vs[1, :]
	M = Matrix{Float64}(LinearAlgebra.I, 4, 4)
	M[1:3, 1:3] = [u1 u2 u3]
	return T * M
end


#//////////////////////////////////////////////////////////////////////////////
function fragment_face_intersection(V::Points, EV::ChainOp, face::Chain; err=LAR_DEFAULT_ERR)
	vs = buildFV(EV, face)     # EV::ChainOp, face::Chain
	retV = Points(undef, 0, 3)
	visited_verts = []
	for i in 1:length(vs)
		o = V[vs[i], :]
		j = i < length(vs) ? i + 1 : 1
		d = V[vs[j], :] - o    # vertex in local coordinates
		if !(-err < d[3] < err)
			alpha = -o[3] / d[3]
			if -err <= alpha <= 1 + err
				p = o + alpha * d
				if -err < alpha < err || 1 - err < alpha < 1 + err
					if !(is_visited_vertex(p, visited_verts))
						push!(visited_verts, p)
						retV = [retV; reshape(p, 1, 3)]
					end
				else
					retV = [retV; reshape(p, 1, 3)]
				end
			end
		end
	end
	vnum = size(retV, 1)
	if vnum == 1
		vnum = 0
		retV = Points(undef, 0, 3)
	end
	enum = (÷)(vnum, 2)
	retEV = spzeros(Int8, enum, vnum)
	for i in 1:enum
		retEV[i, 2*i-1:2*i] = [-1, 1]
	end
	retV, retEV
end



# //////////////////////////////////////////////////////////////////////////////
function fragment_merge_vertices(V::Points, EV::ChainOp, FE::ChainOp; err=LAR_DEFAULT_ERR)
	""" Task to iteratively add new local components to the global 2-skeleton """
	vertsnum = size(V, 1)
	edgenum = size(EV, 1)
	facenum = size(FE, 1)
	newverts = zeros(Int, vertsnum)
	# KDTree constructor needs an explicit array of Float64
	V = Matrix(V)
	kdtree = KDTree(BYROW(V))
	# remove vertices congruent to a single representative
	todelete = []
	i = 1
	for vi in 1:vertsnum
		if !(vi in todelete)
			nearvs = inrange(kdtree, V[vi, :], err)
			newverts[nearvs] .= i
			nearvs = setdiff(nearvs, vi)
			todelete = union(todelete, nearvs)
			i = i + 1
		end
	end
	nV = V[setdiff(collect(1:vertsnum), todelete), :]
	# translate edges to take congruence into account
	edges = Array{Tuple{Int,Int},1}(undef, edgenum)
	oedges = Array{Tuple{Int,Int},1}(undef, edgenum)
	for ei in 1:edgenum
		v1, v2 = EV[ei, :].nzind
		edges[ei] = Tuple{Int,Int}(sort([newverts[v1], newverts[v2]]))
		oedges[ei] = Tuple{Int,Int}(sort([v1, v2]))
	end
	nedges = union(edges)
	# remove edges of zero length
	nedges = filter(t -> t[1] != t[2], nedges)
	nedgenum = length(nedges)
	nEV = spzeros(Int8, nedgenum, size(nV, 1))
	etuple2idx = Dict{Tuple{Int,Int},Int}()
	for ei in 1:nedgenum
		begin
			nEV[ei, collect(nedges[ei])] .= 1
			nEV
		end
		etuple2idx[nedges[ei]] = ei
	end
	for e in 1:nedgenum
		v1, v2 = findnz(nEV[e, :])[1]
		nEV[e, v1] = -1
		nEV[e, v2] = 1
	end
	# compute new faces to take congruence into account
	faces = [[
		map(x -> newverts[x], FE[fi, ei] > 0 ? oedges[ei] : reverse(oedges[ei]))
		for ei in FE[fi, :].nzind
	] for fi in 1:facenum]

	visited = []
	function filter_fn(face)
		verts = []
		map(e -> verts = union(verts, collect(e)), face)
		verts = Set(verts)
		if !(verts in visited)
			push!(visited, verts)
			return true
		end
		return false
	end
	nfaces = filter(filter_fn, faces)
	nfacenum = length(nfaces)
	nFE = spzeros(Int8, nfacenum, size(nEV, 1))
	for fi in 1:nfacenum
		for edge in nfaces[fi]
			ei = etuple2idx[Tuple{Int,Int}(sort(collect(edge)))]
			nFE[fi, ei] = sign(edge[2] - edge[1])
		end
	end
	return Points(nV), nEV, nFE
end

# //////////////////////////////////////////////////////////////////////////////
function fragment_single_face(V, EV::ChainOp, FE::ChainOp, sp_idx, sigma)
	vs_num = size(V, 1)
	# 2D transformation of `sigma` face
	sigmavs = (abs.(FE[sigma:sigma, :])*abs.(EV))[1, :].nzind # sigma vertex indices
	sV = V[sigmavs, :]
	sEV = EV[FE[sigma, :].nzind, sigmavs]
	M = fragment_submanifold_mapping(sV)
	tV = ([V ones(vs_num)]*M)[:, 1:3]  # folle convertire *tutti* i vertici
	sV = tV[sigmavs, :]
	# `sigma` face intersection with faces in `sp_idx[sigma]`, i.e., in `bigpi`
	for i in sp_idx[sigma]
		tmpV, tmpEV = fragment_face_intersection(tV, EV, FE[i, :]) # va in loop qui dentro
		sV, sEV = skel_merge(sV, sEV, tmpV, tmpEV)
	end
	# computation of 2D arrangement of sigma face
	sV = sV[:, 1:2]
	nV, nEV, nFE = planar_arrangement(BYCOL(sV), sEV, sparsevec(ones(Int8, length(sigmavs))))
	nV = BYROW(nV)
	nvsize = size(nV, 1)
	# return each 2D complex in 3D
	nV = [nV zeros(nvsize) ones(nvsize)] * inv(M)[:, 1:3]
	return nV, nEV, nFE
end


# //////////////////////////////////////////////////////////////////////////////
function fragment_compute_copFE(V::Points, copFV::ChainOp, copEV::ChainOp; convex=true::Bool, exterior=false::Bool)::ChainOp

	temp = copFV * copEV'
	I, J, Val = Int64[], Int64[], Int8[]
	for j = 1:size(temp, 2)
		for i = 1:size(temp, 1)
			if temp[i, j] == 2
				push!(I, i)
				push!(J, j)
				push!(Val, 1)
			end
		end
	end
	copFE = SparseArrays.sparse(I, J, Val)
	if !convex
		copFE = fix_redundancy(copFE, copFV, copEV)
	end

	EV = [findnz(copEV[k, :])[1] for k = 1:size(copEV, 1)]
	copEV = sparse(cop_coboundary_0(EV))
	for f = 1:size(copFE, 1)
		chain = findnz(copFE[f, :])[1]#	dense
		cycle = spzeros(Int8, copFE.n)#	sparse

		edge = findnz(copFE[f, :])[1][1]
		sign = 1
		cycle[edge] = sign
		chain = setdiff(chain, edge)
		while chain != []
			boundary = sparse(cycle') * copEV
			_, vs, vals = findnz(dropzeros(boundary))

			rindex = vals[1] == 1 ? vf = vs[1] : vf = vs[2]
			r_boundary = spzeros(Int8, copEV.n)#	sparse
			r_boundary[rindex] = 1
			r_coboundary = copEV * r_boundary
			r_edge = intersect(findnz(r_coboundary)[1], chain)[1]
			r_coboundary = spzeros(Int8, copEV.m)#	sparse
			r_coboundary[r_edge] = EV[r_edge][1] < EV[r_edge][2] ? 1 : -1

			lindex = vals[1] == -1 ? vi = vs[1] : vi = vs[2]
			l_boundary = spzeros(Int8, copEV.n)#	sparse
			l_boundary[lindex] = -1
			l_coboundary = copEV * l_boundary
			l_edge = intersect(findnz(l_coboundary)[1], chain)[1]
			l_coboundary = spzeros(Int8, copEV.m)#	sparse
			l_coboundary[l_edge] = EV[l_edge][1] < EV[l_edge][2] ? -1 : 1

			if r_coboundary != -l_coboundary  # false iff last edge
				# add edge to cycle from both sides
				rsign = rindex == EV[r_edge][1] ? 1 : -1
				lsign = lindex == EV[l_edge][2] ? -1 : 1
				cycle = cycle + rsign * r_coboundary + lsign * l_coboundary
			else
				# add last (odd) edge to cycle
				rsign = rindex == EV[r_edge][1] ? 1 : -1
				cycle = cycle + rsign * r_coboundary
			end
			chain = setdiff(chain, findnz(cycle)[1])
		end
		for e in findnz(copFE[f, :])[1]
			copFE[f, e] = cycle[e]
		end
	end

	pdim=size(V, 1)

	if exterior &&  pdim == 2
		# put matrix in form: first row outer cell; with opposite sign )
		V_row  = BYROW(V)

		outer = get_external_cycle(V_row, copEV, copFE)
		copFE = [-copFE[outer:outer, :]; copFE[1:outer-1, :]; copFE[outer+1:end, :]]
		# induce coherent orientation of matrix rows (see examples/orient2d.jl)
		for k = 1:size(copFE, 2)
			spcolumn = findnz(copFE[:, k])
			if sum(spcolumn[2]) != 0
				row = spcolumn[1][2]
				sign = spcolumn[2][2]
				copFE[row, :] = -sign * copFE[row, :]
			end
		end
		return copFE
	else
		return copFE
	end
end


# ///////////////////////////////////////////////////////////
function fragment_all_faces(lar::Lar)

	V=lar.V
	EV=lar.C[:EV]
	FV=lar.C[:FV]

	# scrgiorgio: if I use here simply `lar2cop` it does not work
	# the cop_XXX function do some magic with orientation
	if true
		copEV = cop_coboundary_0(EV)
		copFE = fragment_compute_copFE(V, lar2cop(FV), lar2cop(EV), convex=true, exterior=false)
	else
		copEV = lar2cop(EV)
		copFV = lar2cop(FV)
		copFE = (copFV * copEV') .÷ Int8(2)
	end
	
	# historically arrangement works internally by using by-row vertices
	V_row = BYROW(V)
	fs_num = size(copFE, 1)
	# strange but necessary cycle of computations to get FV::Cells algebraically
	FV = (abs.(copFE) * abs.(copEV)) .÷ 2
	FV = convert(ChainOp, FV)
	sp_idx = spaceindex(V_row, cop2lar(FV))
	rV = Points(undef, 0, 3)
	rEV = SparseArrays.spzeros(Int8, 0, 0)
	rFE = SparseArrays.spzeros(Int8, 0, 0)
	depot_V = Array{Array{Float64,2},1}(undef, fs_num)
	depot_EV = Array{ChainOp,1}(undef, fs_num)
	depot_FE = Array{ChainOp,1}(undef, fs_num)
	for sigma in 1:fs_num
		print(sigma, "/", fs_num, "\r")
		nV, nEV, nFE = fragment_single_face(V_row, copEV, copFE, sp_idx, sigma)
		depot_V[sigma] = nV
		depot_EV[sigma] = nEV
		depot_FE[sigma] = nFE
	end
	rV = vcat(depot_V...)
	rEV = SparseArrays.blockdiag(depot_EV...)
	rFE = SparseArrays.blockdiag(depot_FE...)
	rV, rcopEV, rcopFE = fragment_merge_vertices(rV, rEV, rFE)
	return rV, rcopEV, rcopFE
end

# ////////////////////////////////////////////////////////////////////////
function lar_tgw(lar::Lar; debug_mode=false)::Cells

  # Newell's method, works for concave too and it is oriented
  function compute_oriented_newell_normal(loop::AbstractPointsNd)::PointNd
    n = [0.0, 0.0, 0.0]
    for (I, (x1, y1, z1)) in enumerate(loop)
        x2, y2, z2 = loop[mod1(I+1, length(loop))]
        n[1] += (y1 - y2) * (z1 + z2)
        n[2] += (z1 - z2) * (x1 + x2)
        n[3] += (x1 - x2) * (y1 + y2)
    end
    return n / LinearAlgebra.norm(n)
  end

  function cycle_normal(V, cycle)
    return compute_oriented_newell_normal([V[:,a] for (a, b) in cycle])
  end

  function cycle_dimension(V, cycle)
    a,b=cycle[1]
    m,M=V[:,a],V[:,a]
    for (a,b) in cycle
      m=[min(m[I],it) for (I,it) in enumerate(V[:,a])]
      M=[max(M[I],it) for (I,it) in enumerate(V[:,a])]
    end
    return prod([M[I]-m[I] for I in 1:3 if (M[I]-m[I])>0]) 
  end

	# tgw angle computation
	function compute_oriented_angle(a::PointNd, b::PointNd, c::PointNd)::Float64
		a = a / LinearAlgebra.norm(a)
		b = b / LinearAlgebra.norm(b)
		angle = atan(dot(cross(a, b), c), dot(a, b))

		ret=rad2deg(angle)

		# I think it will not be a problem if the first face I see is 180 or -180 it means they are two flat faces connected each others
		# @assert(ret!=180.0 || ret==-180.0) # ambiguity otherwise
		return ret
	end

  V=lar.V
  FE=lar.C[:FE]
  EV=lar.C[:EV]
  num_faces=length(FE)

  # compute cycles and normals (note: a normal is unique for a face and shared by all its cycles)
  begin
    cycles,normals=[],[]
    for (A,fe) in enumerate(FE)
      face_cycles=find_vcycles(Cells([EV[E] for E in fe])) 
      if length(face_cycles)>1
        # need to (1) find outer main cycle (2) reverse all other cycles  because they are holes and 
        # need to be traversed in reversed order (even if the normal is the same)
        face_cycles=collect(sort(face_cycles,by=cycle->cycle_dimension(V, cycle),rev=true))
        n1=cycle_normal(V, face_cycles[1])
        for C in 2:length(face_cycles)
          n2=cycle_normal(V, face_cycles[C])
          face_cycles[C]=dot(n1,n2)>=0 ? reverse_cycle(face_cycles[C]) : face_cycles[C]
        end
      end
      push!(cycles, face_cycles)
      push!(normals, cycle_normal(V, face_cycles[1]))
      #println("$(A) NORMAL $(normals[end]) face_cycles=$(face_cycles)")
    end 
  end

  # compute connections
  begin
    connections=Dict{Int,Dict}()
    for (A, Acycles) in enumerate(cycles)
      connections[+A],connections[-A]=Dict(), Dict()
      #println("Face $(A) Acycles=$(Acycles) normal=$(normals[A])")
      for Acycle in Acycles
        for (a, b) in Acycle
          adjacent_pos, adjacent_neg = Vector{Any}(),Vector{Any}()
          for (B, Bcycles) in enumerate(cycles)
            if A==B  continue  end  # same face, skip
            # println("  B=$(B) Bcycles=$(Bcycles)")
            for Bcycle in Bcycles
              if     ([b,a] in Bcycle)  Bsign=+1 # they are correct: in opposite direction 
              elseif ([a,b] in Bcycle)  Bsign=-1 # i must go in the opposite direction
              else  continue  end
              # println("  Found adj cycle edge=$([a,b]) with face $(B) Bsign=$(Bsign)" )
              angle_pos=compute_oriented_angle(+1.0 .* normals[A], +Bsign .* normals[B],   V[:,b] - V[:,a]) # from a->b
              angle_neg=compute_oriented_angle(-1.0 .* normals[A], -Bsign .* normals[B],   V[:,a] - V[:,b]) # from b->a
              push!(adjacent_pos, (face=(+Bsign) * B, angle=angle_pos))
              push!(adjacent_neg, (face=(-Bsign) * B, angle=angle_neg))
            end
          end
          if length(adjacent_pos)==0
            println("PROBLEM WITH edge [$(a),$(b)] only one face incident")
            @assert(false)
          end
          connections[+A][normalize_cell([a,b])]=collect(sort(adjacent_pos,by=value->value.angle))  
          connections[-A][normalize_cell([a,b])]=collect(sort(adjacent_neg,by=value->value.angle))
        end
      end
    end
  end

  # topology checks on connections (if A picks B then B must pick A)
  begin
    for (A,  adjacent_per_edge) in connections
      for (ev, adjacents) in adjacent_per_edge
        for adj in adjacents
          B=adj.face
          @assert B in [it.face for it in connections[A][ev]]
          @assert A in [it.face for it in connections[B][ev]]
        end
      end
    end
  end

  # find best connections by angles (TGW)
  begin
    best_connections=Dict{Int, Cell}()
    for (A,  adjacent_per_edge) in connections
      best_connections[A]=remove_duplicates([adjacents[1].face for (ev, adjacents) in adjacent_per_edge])
    end
  end

  #for (k,v) in best_connections
  #  println(k," ",v)
  #end

  @assert(all([abs(F)>=1 && abs(F)<=num_faces for F in keys(best_connections)]))
  atoms=lar_connected_components(collect(keys(best_connections)), cur -> best_connections[cur])
  # @show(atoms)

  # atoms topology checks
  begin
    for F in keys(connections)
      # one signed face should be in only one atom
      @assert(length([atom for atom in atoms if F in atom])==1)
      # if there is +F in one atom, there cannot be -F
      @assert(length([atom for atom in atoms if +F in atom && -F in atom])==0)
    end
  end

  # I do not need the sign anymore (could be useful for debugging)
  atoms=[collect(sort([abs(jt) for jt in it])) for it in atoms]

  # @show(atoms)

  # topology check
  begin
    num_full_cell_per_face= Dict{Int,Int}()
    for (A, atom) in enumerate(atoms)

      # each atom must not contain the same face twice
      if length(Set(atom))!=length(atom)
        @assert(false)
      end

      for F in atom
        if !haskey(num_full_cell_per_face,F) num_full_cell_per_face[F]=0 end
        num_full_cell_per_face[F]+=1
      end
    end

    # @show(num_full_cell_per_face)

    # no hanging face, and a face cannot be shared by more than 2-full-dim cells
    # note: consider the outer atom too
    @assert(all([v==2 for (k,v) in num_full_cell_per_face]))

  end

  return atoms
end

# //////////////////////////////////////////////////////////////////////////////
function ARRANGE3D(lar::Lar; debug_mode=false)::Lar
	""" Main function of arrangement pipeline """
	rV, rcopEV, rcopFE=fragment_all_faces(lar)

  EV = cop2lar(rcopEV)
  FE = cop2lar(rcopFE) 
  FV = [union(CAT([EV[E] for E in fe])) for fe in FE]

  ret = Lar(BYCOL(rV),Dict(
		:EV => EV, 
		:FE => FE, 
		:FV => FV
	))

	if debug_mode
		show_debug(ret)
		VIEWCOMPLEX(ret, explode=[1.2,1.2,1.2], show=["V","FV","Vtext"], title="arrange3d / 3d ALL faces")
	end  

	ret=SIMPLIFY(ret) # not sure this is needed

	ret.C[:FV]=compute_FV(ret)
	ret.C[:CF]=lar_tgw(ret, debug_mode=debug_mode)

	# ret.C[:CV]=compute_CV(ret) dont think this is needed
	CHECK(ret)

	if debug_mode
		@show(ret)
		VIEWCOMPLEX(ret, show=["CV"], explode=[1.2,1.2,1.2], title="arrange3d / ALL atoms")
	end

	ret=SIMPLIFY(ret)
	# @show(ret)
	return ret
end
export ARRANGE3D







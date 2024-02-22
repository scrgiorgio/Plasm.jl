function vertCongruenceAA(W)
	err, i, todelete, vclasses = 10^-6, 1, [], []
	verts = convert(Lar.Points, W')
	kdtree = NearestNeighbors.KDTree(verts)
	newverts = zeros(Int, size(verts,2))
	for vi in 1:size(verts,2)
		if !(vi in todelete)
			nearvs = NearestNeighbors.inrange(kdtree, verts[:,vi], err)
			push!(vclasses,nearvs)
			newverts[nearvs] .= i
			nearvs = setdiff(nearvs, vi)
			todelete = union(todelete, nearvs)
			i += 1
		end
	end
	V = zeros(3,length(vclasses))
	for (k,class) in enumerate(vclasses)
			V[:,k] = sum(W[class,:],dims=1)/length(class)
	end
	return V, vclasses
end


function cellCongruenceAA(Delta,inclasses)
	cellarray = Lar.cop2lar(Delta)
	new_e = Array{Int64,1}(undef,size(Delta,2))
	for (k,class) in enumerate(inclasses)
		for e in class
			new_e[e] = k
		end
	end
  cells = [map(e->new_e[e], face) for face in cellarray]
  outclasses = DefaultOrderedDict{Array{Int64,1},Array{Int64,1}}([])
  for (k,face) in enumerate(cells)
    if outclasses[face] == []
    	outclasses[face] = [k]
    else
    	append!(outclasses[face],[k])
    end
  end
	FEnew = sort(collect(keys(outclasses)))
  outclasses = sort(collect(values(outclasses)))
  return FEnew,outclasses
end


function chainCongruenceAA(W, T)
	V, vclasses = vertCongruenceAA(W)
	EV, eclasses = cellCongruenceAA(T[1],vclasses)
	FE, fclasses = cellCongruenceAA(T[2],eclasses)
	#@show size(V,2) - size(EV,1) + size(FE,1)
	return V,EV,FE
end
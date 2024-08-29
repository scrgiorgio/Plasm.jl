using Plasm,SparseArrays,NearestNeighbors

export cop2lar, coboundary_1, vcongruence, cellcongruence, chaincongruence


#///////////////////////////////////////////////////////////////////////////////
function coboundary_1(FV::Array{Array{Int64,1},1}, EV::Array{Array{Int64,1},1}) # (::Cells, ::Cells)
copFV =  lar2cop(FV)
I,J,Val = findnz( lar2cop(EV))
copVE = sparse(J,I,Val)
triples = hcat([[i,j,1]  for (i,j,v)  in zip(findnz(copFV * copVE)...) if v==2]...)
I,J,Val = triples[1,:], triples[2,:], triples[3,:]
Val = convert(Array{Int8,1},Val)
copFE = sparse(I,J,Val)
return copFE
end

#///////////////////////////////////////////////////////////////////////////////
function chaincongruence(W, Delta_0, Delta_1; epsilon=1e-6)
	V, vclasses = vcongruence(W, epsilon); map(sort!,vclasses)
	EV = cellcongruence(Delta_0, vclasses, dim=1)
	EV = sort(union(map(sort!, EV)))
	FV = cellcongruence(Delta_1 * Delta_0, vclasses, dim=2)
	FV = sort(union(map(sort!, FV)))
	copFE = coboundary_1(FV:: Cells, EV:: Cells)
	FE = cop2lar(copFE)
	return V, EV, FV, FE
end

#///////////////////////////////////////////////////////////////////////////////
function vcongruence(V::Matrix,epsilon)
    vclasses, visited = [], []
    kdtree = NearestNeighbors.KDTree(V);
    for vidx = 1 : size(V, 2) if !(vidx in visited)
       nearvs = NearestNeighbors.inrange(kdtree, V[:,vidx], epsilon)
       push!(vclasses, nearvs)
			 append!(visited, nearvs) end
    end
    W = hcat([sum(V[:,class], dims=2)/length(class) for class in vclasses]...)
    return W, vclasses
end


#///////////////////////////////////////////////////////////////////////////////
function cellcongruence(Delta, vclasses; dim)  # ==> OK !!!
  # cells  with old vertices
	cellarray =  cop2lar(Delta) 
	# conversion array old -> new vertex
	newvert = Vector(undef, size(Delta,2))
  for (k, class) in enumerate(vclasses) 
  	for v in class 
  	  newvert[v] = k
  	end
  end	
	# conversione delle celle da vecchi a nuovi indici di vertice
	newcells = Vector{Int64}[]
	for cell in cellarray
		newcell = Int64[]
		for v in cell
			push!(newcell, newvert[v])
		end
		length(newcell) > dim ? push!(newcells, newcell) : break
	end
	# remove duplicate cells
	newcells = union(map(sort,newcells))
		# remove duplicate vertices
		newcells = map(union,newcells)
	# remove empty cells
	outcells = filter(x -> !(length(x) <= dim || all(y->y==x[1], x)), newcells)
end


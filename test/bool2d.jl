using LinearAlgebraicRepresentation, SparseArrays
Lar = LinearAlgebraicRepresentation;
using IntervalTrees,LinearAlgebra
# using Revise, OhMyREPL

using DataStructures, Base

function LAR(obj::Hpc)::Lar
	 obj=ToGeometry(obj)
   V  = obj.points;
   EV = obj.edges;
   FV = obj.facets;
   V,FV,EV = simplifyCells(hcat(V...),FV) # !!!!  simplifyCells(hcat(V...),FV,EV);
   if !(FV == [])
      FV = union(FV)
      FF = CSC(FV) * CSC(FV)'
      edges = filter(x->x[1]<x[2] && FF[x...]==2, collect(zip(findnz(FF)[1:2]...)))
      EW = sort!(collect(Set([FV[i] ∩ FV[j] for (i,j) in edges])))
   elseif all(length(x)==2 for x in EV) 
      return Plasm.Lar(1,V, Dict(:EV=>EV))
   end
   return Plasm.Lar(V, Dict(:FV=>FV, :EV=>EW))
end


function simplifyCells(V,CV)
	PRECISION = 14
	vertDict = DefaultOrderedDict{Vector{Float64}, Int64}(0)
	index = 0
	W = Vector{Float64}[]
	FW = Vector{Int64}[]

	for incell in CV
		outcell = Int64[]
		for v in incell
			vert = V[:,v]
			key = map(Base.truncate(PRECISION), vert)
			if vertDict[key]==0
				index += 1
				vertDict[key] = index
				push!(outcell, index)
				push!(W,key)
			else
				push!(outcell, vertDict[key])
			end
		end
		append!(FW, [[Set(outcell)...]])
	end
	return hcat(W...), filter(x->!(LEN(x)<2), FW)
end

function simplifyCells(V,FV,EV)
   V,FV = simplifyCells(V,CV)
   V,EV = simplifyCells(V,EV)
   return V,FV,EV
end


# //////////////////////////////////////////////////////////////////////////////
function settestpoints2d(W,EV,FV,f, copEV,copFE) # W by rows
	e = findnz(copFE[f,:])[1][1] # first edge of face f
	v1,v2 = findnz(copEV[e,:])[1] # two (global) verts incident on in
    t = W[v2,:] - W[v1,:]
    n = [-t[2],t[1]]
	p0 = (W[v1,:] + W[v2,:]) ./ 2
	ϵ = 1.0e-4
	ptest1 = p0 + ϵ * n
	ptest2 = p0 - ϵ * n
	return ptest1, ptest2
end

# //////////////////////////////////////////////////////////////////////////////
function getinternalpoint2d(W,EV,FV, f, copEV,copFE) # W by rows
	#edges for v1=FV[1][1]
	ptest1, ptest2 = settestpoints2d(W,EV,FV, f, copEV,copFE)
    edges = [findnz(copEV[e,:])[1] for e in findnz(copFE[f,:])[1]]
    V = convert(Points,W')
    classify = pointInPolygonClassification(V,edges)
    if classify(ptest1) == "p_in"
        return ptest1
    elseif classify(ptest2) == "p_in"
        return ptest2
    else
        error("classifying inner point in face $f")
    end
end

# //////////////////////////////////////////////////////////////////////////////
function chainbasis2polygons(V,copEV,copFE)
	FE = [findnz(copFE[k,:])[1] for k=1:copFE.m]
	EV = [findnz(copEV[k,:])[1] for k=1:copEV.m]

	FEs = Array{Int64,1}[]
	EVs = Array{Array{Int64,1},1}[]
	FVs = Array{Int64,1}[]
	for f=1:copFE.m
		push!( FEs, collect(Set(cat([e for e in FE[f]]))) )
		# edges in EVs are aggregated by face, in order to answer point-classifications
		push!( EVs, [EV[e] for e in FE[f]]  )
		push!( FVs, collect(Set(cat([EV[e] for e in FE[f]])))  )
	end
	polygons = collect(zip(EVs,FVs,FEs)) # (EV,FV,FFE)s for polygon
	W = convert(Points,V')
	return W,polygons,FE
end


# //////////////////////////////////////////////////////////////////////////////
function internalpoints2d(W,copEV,copFE) # W by rows
	# transform each 2-cell in a solid (via Lar model)
   #----------------------------------------------------------------------------
	U,pols,FE = chainbasis2polygons(W,copEV,copFE) # V by rows
	# compute, for each `pol` (3-cell) in `pols`, one `internalpoint`.
	#----------------------------------------------------------------------------
	internalpoints = []
	for f=1:length(pols)
		(EV,FV,FE) = pols[f]
		#GL.VIEW([ GL.GLFrame, GL.GLLines(V,EV) ]);
		internalpoint = getinternalpoint2d(W,EV,FV, f, copEV,copFE)
		push!(internalpoints,internalpoint)
	end
	return internalpoints
end

"""
	testinternalpoint2d(V::Points, EV::Cells, FV::Cells)

"""
function testinternalpoint2d(listOfModels)
	function testinternalpoint0(testpoint)
		intersectedfaces = Int64[]
		depot = []
		# actual containment test of ray point in faces within depot
		for (k,model) in enumerate(listOfModels)
            verts,edges = model
			classify = pointInPolygonClassification(verts,edges)
			inOut = classify(testpoint)
			if inOut == "p_in"
				push!(intersectedfaces,k)
			end
		end
		return intersectedfaces
	end
	return testinternalpoint0
end


################################################################################

function bool2d(assembly::Hpc)
	# input of affine assembly
	# input of affine assembly
   ensemble = LAR(assembly)
   V,FV,EV = new2old(ensemble)
	#----------------------------------------------------------------------------
	# V,EV = Lar.struct2lar(assembly) #TODO proper different method
	V,EV = Lar.struct2lar(assembly)
	cop_EW = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells));
	W = convert(Lar.Points, V');
	# generate the 3D space arrangement
	#----------------------------------------------------------------------------
	W, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)
	#V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF)
	innerpoints = Lar.internalpoints2d(W,copEV,copFE[1:end,:])
	# associate internal points to 2-cells
	#----------------------------------------------------------------------------
	listOfModels = Lar.evalStruct(assembly)
	inputfacenumbers = [length(listOfModels[k][2]) for k=1:length(listOfModels)]
	# test input data for containment of reference points
	#----------------------------------------------------------------------------
	boolmatrix = BitArray(undef, length(innerpoints), length(listOfModels))
    containmenttest = Lar.testinternalpoint2d(listOfModels)
	for (k,point) in enumerate(innerpoints) # k runs on columns
		cells = containmenttest(point) # contents of columns
		#println(k," ",faces)
		for l in cells
			boolmatrix[k,l] = 1
		end
	end
	return W, copEV, copFE, boolmatrix
end


"""
	boolops( assembly::Hpc, op::Symbol )

User interface, where `op` symbol ``in`` {`:|`, `:&`, `:-`}. Return a pair V,EV
# Example
```julia
assembly = STRUCT(wall, openings)
V,EV = boolops(assembly, :-)
```
"""
# //////////////////////////////////////////////////////////////////////////////
function boolops(assembly::Lar.Struct, op::Symbol)
    W, copEV, copFE, boolmatrix = Lar.bool2d(assembly)
    boolvars = [boolmatrix[:,k] for k=1:size(boolmatrix,2)]
    if eval(op) == -
        solution = boolvars[1] .& .!(.|(boolvars[2:end]...))
    else
        operator = eval(op)
        solution = operator.(boolvars...)
    end
    chain1d = sparse(copFE') * Int.(solution)
    EV = Lar.cop2lar(copEV)
    EVop = [ev for (k,ev) in enumerate(EV) if abs(chain1d[k])==1 ]
    V = convert(Lar.Points,W')
    return V,EVop
end




################################################################################

using Plasm

# //////////////////////////////////////////////////////////////////////////////
# 2D Boolean example generation (see CAD23 paper)
n,m = 1,1
square = TYPE(Hpc(CUBOIDGRID([n,m])), "solid")

assembly = STRUCT( 
    STRUCT( T(1,2)(-√2/4, -√2/2 ), R(1,2)(π/4), square ),
    STRUCT( T(1,2)( √2/4, -√2/2 ), R(1,2)(π/4), square ))
    
VIEW(assembly)


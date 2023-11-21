using Plasm

const Points = Matrix
const Cells = Vector{Vector{Int}}
#const Cell = SparseVector{Int8, Int}
#const Chain = SparseVector{Int8,Int}
#const ChainOp = SparseMatrixCSC{Int8,Int}


# //////////////////////////////////////////////////////////////////////////////

# method with Hpc data

function EXPLODER(V,FVs,sx=1.5,sy=1.5,sz=1.5)
	outVerts, outCells = [],[]
	pols = Hpc[]
	for FV in FVs
		#vertidx = sort(collect(Set(cat(FV,dims=1))))
		vertidx = if !(FV==Any[]) sort(union(FV...))
		    else break end
		vcell = V[:,vertidx]
		vdict = Dict(zip(vertidx,1:length(vertidx)))

		center = sum(vcell,dims=2)/size(vcell,2)
		scaled_center = size(center,1)==2 ? center .* [sx;sy] : center .* [sx;sy;sz]
		translation_vector = scaled_center - center
		cellverts = vcell .+ translation_vector
		newcells = [[vdict[v] for v in cell] for cell in FV]
      U =  [cellverts[:,k] for k=1:size(cellverts,2)]
		pol = MKPOL( U, convert(Cells,newcells) )
		push!(pols, pol)
	end
	return STRUCT(pols)
end

# method with Hpc data
function EXPLODER(pols::Vector{Hpc},sx=1.5,sy=1.5,sz=1.5)::Vector{Hpc}
   n = LEN(pols)
   #vect = [(-1 .+[sx,sy,sz]*5) for k=1:n]
   vect = [(-1 .+[sx,sy,sz]*5) for k=1:n]
   newpols = [T([1,2,3])(vect[k])(cubes[k]) for k=1:n]
end

# //////////////////////////////////////////////////////////////////////////////
# Functio to color exploded assemblies ///////////// TODO //////////////////////

#function GLExplode(V,FVs,sx=1.2,sy=1.2,sz=1.2,colors=1,alpha=0.2::Float64)
#	assembly = GL.explodecells(V,FVs,sx,sy,sz)
#	meshes = Any[]
#	for k=1:length(assembly)
#		if assembly[k] ≠  Any[]
#			# Lar model with constant lemgth of cells, i.e a GRID object !!
#			V,FV = assembly[k]
#			col = GL.Point4d(1,1,1,1)
#			# cyclic color + random color components
#			if colors == 1
#				col = GL.COLORS[1]
#			elseif 2 <= colors <= 12
#				col = GL.COLORS[colors]
#			else # colors > 12: cyclic colors w random component
#				col = GL.COLORS[(k-1)%12+1] - (rand(Float64,4)*0.1)
#			end
#			#col *= alpha
#			push!(meshes, GL.GLGrid(V,FV,col,alpha) )
#		end
#	end
#	return meshes
#end

# //////////////////////////////////////////////////////////////////////////////


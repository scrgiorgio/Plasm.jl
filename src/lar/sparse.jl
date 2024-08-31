

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
function cop_boundary_0(EV::Cells)::ChainOp
	copEV = lar2cop(EV)
	copVE = copEV'
	for (E,ev) in enumerate(EV)
		v1,v2=ev
		copVE[v1, E] = -1 # from +1 -> -1
	end
	copEV=LinearAlgebra.transpose(copVE)
	return convert(ChainOp,copEV)
end
export cop_boundary_0

# //////////////////////////////////////////////////////////////////////////////
"""From (EV,FE) to EV"""
function FV2EV(copEV::ChainOp, copFE::ChainOp)
	EV = cop2lar(copEV) 
	FE = cop2lar(copFE)
	return [[EV[e] for e in fe] for fe in FE]
end
export FV2EV

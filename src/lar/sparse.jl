

const Chain        = SparseVector{Int8,Int}
const ChainOp      = SparseMatrixCSC{Int8,Int}
const ChainComplex = Vector{ChainOp}
export Chain, ChainOp, ChainComplex

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


""" converte sparse represantation to dense"""
function cop2lar(cop::ChainOp)::Cells
	[findnz(cop[k, :])[1] for k = 1:size(cop, 1)]
end
export cop2lar

function boundary_1(EV::Cells)::ChainOp
	out = lar2cop(EV)'
	for e = 1:length(EV)
		out[EV[e][1], e] = -1
	end
	return out
end
export boundary_1

function coboundary_0(EV::Cells)::ChainOp
	return convert(ChainOp, LinearAlgebra.transpose(boundary_1(EV::Cells)))
end
export coboundary_0

# //////////////////////////////////////////////////////////////////////////////
"""Find EV from EV FE"""
function FV2EVs(copEV::ChainOp, copFE::ChainOp)
	EV = [findnz(copEV[k, :])[1] for k = 1:size(copEV, 1)]
	FE = [findnz(copFE[k, :])[1] for k = 1:size(copFE, 1)]
	EVs = [[EV[e] for e in fe] for fe in FE]
	return EVs
end
export FV2EVs

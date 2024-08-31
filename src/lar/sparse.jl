

const Chain        = SparseVector{Int8,Int}
const ChainOp      = SparseMatrixCSC{Int8,Int}
const ChainComplex = Vector{ChainOp}
export Chain, ChainOp, ChainComplex

function characteristicMatrix(FV::Cells)::ChainOp
	I, J, V = Int64[], Int64[], Int8[]
	for f = 1:length(FV)
		for k in FV[f]
			push!(I, f)
			push!(J, k)
			push!(V, 1)
		end
	end
	return sparse(I, J, V)
end
export characteristicMatrix

""" converte dense represantation to sparse"""
function lar2cop(CV::Cells)::ChainOp
	I = Int[]
	J = Int[]
	Value = Int8[]
	for k = 1:size(CV, 1)
		n = length(CV[k])
		append!(I, k * ones(Int, n))
		append!(J, CV[k])
		append!(Value, ones(Int, n))
	end
	return SparseArrays.sparse(I, J, Value)
end
export lar2cop

""" converte sparse represantation to dense"""
function cop2lar(cop::ChainOp)::Cells
	[findnz(cop[k, :])[1] for k = 1:size(cop, 1)]
end
export cop2lar

function boundary_1(EV::Cells)::ChainOp
	out = characteristicMatrix(EV)'
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

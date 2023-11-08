using SparseArrays
#= Sparse matrix generation from a set of {triplets} 
(i,j,x) of non-zero, given as a triple (I,J,X) of arrays =#
"""
Characteristic sparse binary matrix of CV collection of h-cells.
"""
function K( CV ) # CV => Cells defined by their Vertices
	I = vcat( [ [k for h in CV[k]] for k=1:length(CV) ]...)		
	# vcat maps arrayofarrays to single array
	J = vcat( CV...)	 
	# splat transforms the CV array elements to vcat arguments
	X = Int8[1 for k=1:length(I)] 	
	# Type Int8 defines the memory map of array elements
	return SparseArrays.sparse(I,J,X)	
end
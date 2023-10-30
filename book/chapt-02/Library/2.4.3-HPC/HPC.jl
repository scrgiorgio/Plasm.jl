mutable struct HPC
    d::Int # intrinsic dimension
    m::Int # embedding dimension (rows of V)
    n::Int # number of vertices  (columns of V)
    V::Matrix{Float64} # object geometry
    C::Dict{Symbol, AbstractArray} # object topology (C for cells)
    # inner constructors
    HPC() = new( -1, 0, 0, Matrix{Float64}(undef,0,0), Dict{Symbol, AbstractArray}() )
    HPC(m::Int,n::Int) = new( m,m,n, Matrix(undef,m,n), Dict{Symbol,AbstractArray}() )
    HPC(d::Int,m::Int,n::Int) = new( d,m,n, Matrix(undef,m,n), Dict{Symbol,AbstractArray}() ) 
    HPC(V::Matrix) = begin m, n = size(V); 
       new( m,m,n, V, Dict{Symbol,AbstractArray}() ) end
    HPC(V::Matrix,C::Dict) = begin m,n = size(V); new( m,m,n, V, C )  end
    HPC(d::Int,V::Matrix,C::Dict) = begin m,n = size(V); new( d,m,n, V, C )  end
    HPC(d,m,n, V,C) = new( d,m,n, V,C )
end

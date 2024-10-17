
""" Module for integration of polynomials over 3D volumes and surfaces """

using Plasm, LinearAlgebra

export SURFACE, VOLUME, FIRSTMOMENT, SECONDMOMENT, INERTIAPRODUCT, CENTROID, INERTIAMOMENT, M, TTT, II, III

#///////////////////////////////////////////////////////////////
function M(alpha::Int, beta::Int)::Float64
    a = 0
    for l=0:(alpha + 1)
        a += binomial(alpha+1,l) * (-1)^l/(l+beta+1)
    end
    return a/(alpha + 1)
end

#///////////////////////////////////////////////////////////////
""" The main integration routine """
function TTT(
tau::Array{Float64,2}, 
alpha::Int, beta::Int, gamma::Int, 
signedInt::Bool=false)
	vo,va,vb = tau[:,1],tau[:,2],tau[:,3]
	a = va - vo
	b = vb - vo
	s1 = 0.0
	for h=0:alpha
		for k=0:beta
			for m=0:gamma
				s2 = 0.0
				for i=0:h 
					s3 = 0.0
					for j=0:k
						s4 = 0.0
						for l=0:m
							s4 += binomial(m,l) * a[3]^(m-l) * b[3]^l * M( 
								h+k+m-i-j-l, i+j+l )
						end
						s3 += binomial(k,j) * a[2]^(k-j) * b[2]^j * s4
					end
					s2 += binomial(h,i) * a[1]^(h-i) * b[1]^i * s3;
				end
				s1 += binomial(alpha,h) * binomial(beta,k) * binomial(gamma,m) * 			
						vo[1]^(alpha-h) * vo[2]^(beta-k) * vo[3]^(gamma-m) * s2
			end
		end
	end
	c = LinearAlgebra.cross(a,b)
	if signedInt == true
		return s1 * LinearAlgebra.norm(c) * sign(c[3])
	else
		return s1 * LinearAlgebra.norm(c)
	end	
end

#///////////////////////////////////////////////////////////////
""" 
	II(P, alpha::Int, beta::Int, gamma::Int, signedInt=false)

Basic integration function on 2D plane.

# Example  unit 3D triangle
```julia
julia> V = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
3×3 Array{Float64,2}:
 0.0  1.0  0.0
 0.0  0.0  1.0
 0.0  0.0  0.0

julia> FV = [[1,2,3]]
1-element Array{Array{Int64,1},1}:
 [1, 2, 3]

julia> P = V,FV
([0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0], Array{Int64,1}[[1, 2, 3]])

julia> II(P, 0,0,0)
0.5
```
"""
function II(
P, 
alpha::Int, beta::Int, gamma::Int, 
signedInt=false)::Float64
    w = 0
    V, FV = P
    if typeof(FV) == Array{Int64,2}
    	FV = [FV[:,k] for k=1:size(FV,2)]
    end
    for i=1:length(FV)
        tau = hcat([V[:,v] for v in FV[i]]...)
        if size(tau,2) == 3
        	term = TTT(tau, alpha, beta, gamma, signedInt)
        	if signedInt
        		w += term
        	else
        		w += abs(term)
        	end
        elseif size(tau,2) > 3
        	println("ERROR: FV[$(i)] is not a triangle")
        else
        	println("ERROR: FV[$(i)] is degenerate")
        end
    end    
    return w
end

#///////////////////////////////////////////////////////////////
""" 
	III(P, alpha::Int, beta::Int, gamma::Int)::Float64

Basic integration function on 3D space.

# Example # unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
3×4 Array{Float64,2}:
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0

julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]]
4-element Array{Array{Int64,1},1}:
 [1, 2, 4]
 [1, 3, 2]
 [4, 3, 1]
 [2, 3, 4]

julia> P = V,FV
([0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0], 
Array{Int64,1}[[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]])

julia> III(P, 0,0,0)
0.16666666666666674
```
"""
function III(P, alpha::Int, beta::Int, gamma::Int)::Float64
    w = 0
    V, FV = P
    for i=1:length(FV)
        tau = hcat([V[:,v] for v in FV[i]]...)
        vo,va,vb = tau[:,1],tau[:,2],tau[:,3]
        a = va - vo
        b = vb - vo
        c = LinearAlgebra.cross(a,b)
        w += c[1]/LinearAlgebra.norm(c) * TTT(tau, alpha+1, beta, gamma)
    end
    return w/(alpha + 1)
end




#///////////////////////////////////////////////////////////////
"""
	SURFACE(P, signedInt::Bool=false)::Float64

`SURFACE` integral on polyhedron `P`.

# Example # unit cube
```julia
julia> V = [0.0 1.0 0.0 0.0 1.0 1.0 0.0 1.0; 0.0 0.0 1.0 0.0 1.0 0.0 1.0 1.0; 0.0 0.0 0.0 1.0 0.0 1.0 1.0 1.0]

julia> FV = [[1,2,6],[6,4,1],[1,3,5],[5,2,1],[4,7,3],[3,1,4],[2,5,8], [8,6,2],[3,7,8],[8,5,3],[4,6,8],[8,7,4]]

julia> P = V,FV
([0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0], [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]])

julia> SURFACE(P)
6,0
```
"""
function SURFACE(P, signedInt::Bool=false)::Float64
    return II(P, 0, 0, 0, signedInt)
end



#///////////////////////////////////////////////////////////////
"""
	VOLUME(P)::Float64

`VOLUME` integral on polyhedron `P`.

# Example 1 # standard unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
julia> P = V,FV
([0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0], [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]])
julia> VOLUME(P)
0.16666666666666666
```
# Example 1 # standard unit cube
```julia
V = [0.0 1.0 0.0 0.0 1.0 1.0 0.0 1.0; 0.0 0.0 1.0 0.0 1.0 0.0 1.0 1.0; 0.0 0.0 0.0 1.0 0.0 1.0 1.0 1.0]
FV= [[1,2,6],[6,4,1],[1,3,5],[5,2,1],[4,7,3],[3,1,4],[2,5,8],
[8,6,2],[3,7,8],[8,5,3],[4,6,8],[8,7,4]]
cube = V,FV
VOLUME(cube) # => 1.0```
"""
function VOLUME(P)::Float64
    return III(P, 0, 0, 0)
end

#v,(vv,ev,fv,cv) = p.cuboidGrid((1,1,1),true)
#V = hcat([Array{Float64,1}(v[k,:]) for k=1:size(v,1)]...)
#FV = hcat([Array{Int64,1}(fv[k,:]+1) for k=1:size(fv,1)]...)
#EV = hcat([Array{Int64,1}(ev[k,:]+1) for k=1:size(ev,1)]...)
#model1 = Any[V,FV,EV]
#P = V,[FV[:,k] for k=1:size(FV,2)]
#SURFACE(P,false)
#
#
#""" Surface integration """
#function surfIntegration(larModel)
#    V,FV,EV = model
#    FE = crossRelation(FV,EV)
#    if typeof(FV) == Array{Int64,2}
#    	FV = [FV[:,k] for k=1:size(FV,2)]
#    end
#    if typeof(V) == Array{Int64,2}
#    	if size(V,1) == 2
#    		V = vcat(V,zeros(1,size(V,2)))
#    	end
#    end
#    cochain = []
#    triangles = []
#    faceVertPairs = []
#	for face=1:length(FE)
#		push!(faceVertPairs, hcat([EV[:,e] for e in FE[face]]...))
#		row = [faceVertPairs[face][1] for k=1:length(FE[face])]
#		push!(triangles, vcat(faceVertPairs[face],row'))
#        P = V,triangles[face]
#        area = SURFACE(P,false) 
#        push!(cochain,area)
#    end
#    return cochain
#end
    
#    TODO: after having implemented ∂_3/∂_2
#def signedIntSurfIntegration(model,signedInt=False)
#    V,FV,EV = model
#    V = [v+[0.0] if len(v)==2 else v for v in V]
#    cochain = []
#    for face in FV
#        triangles = AA(C(AL)(face[0]))(TRANS([face[1-1],face[2]]))
#        P = V,triangles
#        area = Surface(P,signedInt) 
#        cochain += [area]
#    return cochain
#


#///////////////////////////////////////////////////////////////
""" 
	FIRSTMOMENT(P)::Array{Float64,1}

First moments as terms of the Euler tensor. Remember that the integration algorithm is a boundary integration. Hence the model must be a boundary model. In this case, a 2-complex of triangles. 

# Example # unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
julia> P = V,FV;
julia> FIRSTMOMENT(P)
3-element PointNd:
 0.04166666666666663
 0.041666666666666664
 0.041666666666666685
```
# Example 2 # standard unit cube
```julia
V = [0.0 1.0 0.0 0.0 1.0 1.0 0.0 1.0; 0.0 0.0 1.0 0.0 1.0 0.0 1.0 1.0; 0.0 0.0 0.0 1.0 0.0 1.0 1.0 1.0]
FV= [[1,2,6],[6,4,1],[1,3,5],[5,2,1],[4,7,3],[3,1,4],[2,5,8],
[8,6,2],[3,7,8],[8,5,3],[4,6,8],[8,7,4]]
cube = V,FV
FIRSTMOMENT(cube)
3-element PointNd:
 0.5
 0.5
 0.5
```
"""
function FIRSTMOMENT(P)::Array{Float64,1}
    out = zeros(3)
    out[1] = III(P, 1, 0, 0)
    out[2] = III(P, 0, 1, 0)
    out[3] = III(P, 0, 0, 1)
    return out
end



#///////////////////////////////////////////////////////////////
""" 
	SECONDMOMENT(P)::Array{Float64,1}

Second moments as terms of the Euler tensor.

# Example # unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];

julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];

julia> P = V,FV;

julia> SECONDMOMENT(P)
3-element PointNd:
 0.01666666666666668
 0.016666666666666663
 0.016666666666666663
```
# Example 2 # standard unit cube
```julia
V = [0.0 1.0 0.0 0.0 1.0 1.0 0.0 1.0; 0.0 0.0 1.0 0.0 1.0 0.0 1.0 1.0; 0.0 0.0 0.0 1.0 0.0 1.0 1.0 1.0]
FV= [[1,2,6],[6,4,1],[1,3,5],[5,2,1],[4,7,3],[3,1,4],[2,5,8],
[8,6,2],[3,7,8],[8,5,3],[4,6,8],[8,7,4]]
cube = V,FV
SECONDMOMENT(cube)
3-element PointNd:
 0.3333333333333333
 0.33333333333333326
 0.33333333333333326
```
"""
function SECONDMOMENT(P)::Array{Float64,1}
    out = zeros(3)
    out[1] = III(P, 2, 0, 0)
    out[2] = III(P, 0, 2, 0)
    out[3] = III(P, 0, 0, 2)
    return out
end



#///////////////////////////////////////////////////////////////
""" 
	INERTIAPRODUCT(P)::Array{Float64,1}

Inertia products as terms of the Euler tensor.

# Example # unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
julia> P = V,FV;
julia> INERTIAPRODUCT(P)
3-element PointNd:
 0.008333333333333361
 0.008333333333333331
 0.008333333333333325
```
# Example 2 # standard unit cube
```julia
V = [0.0 1.0 0.0 0.0 1.0 1.0 0.0 1.0; 0.0 0.0 1.0 0.0 1.0 0.0 1.0 1.0; 0.0 0.0 0.0 1.0 0.0 1.0 1.0 1.0]
FV= [[1,2,6],[6,4,1],[1,3,5],[5,2,1],[4,7,3],[3,1,4],[2,5,8],
[8,6,2],[3,7,8],[8,5,3],[4,6,8],[8,7,4]]
cube = V,FV
INERTIAPRODUCT(cube)
3-element PointNd:
 0.24999999999999994
 0.25
 0.25
```
"""
function INERTIAPRODUCT(P)::Array{Float64,1}
    out = zeros(3)
    out[1] = III(P, 0, 1, 1)
    out[2] = III(P, 1, 0, 1)
    out[3] = III(P, 1, 1, 0)
    return out
end



#///////////////////////////////////////////////////////////////
""" 
	CENTROID(P)::Array{Float64,1}

Barycenter or `CENTROID` of polyhedron `P`.

# Example # unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
julia> P = V,FV;
julia> CENTROID(P)
3-element PointNd:
 0.24999999999999978
 0.25
 0.2500000000000001
```
# Example 2 # standard unit cube
```julia
V = [0.0 1.0 0.0 0.0 1.0 1.0 0.0 1.0; 0.0 0.0 1.0 0.0 1.0 0.0 1.0 1.0; 0.0 0.0 0.0 1.0 0.0 1.0 1.0 1.0]
FV= [[1,2,6],[6,4,1],[1,3,5],[5,2,1],[4,7,3],[3,1,4],[2,5,8],
[8,6,2],[3,7,8],[8,5,3],[4,6,8],[8,7,4]]
cube = V,FV
CENTROID(cube)
3-element PointNd:
 0.5
 0.5
 0.5
```
"""
function CENTROID(P)::Array{Float64,1}
	return FIRSTMOMENT(P)./VOLUME(P)
end



#///////////////////////////////////////////////////////////////
""" 
	INERTIAMOMENT(P)::Array{Float64,1}

Inertia moments  of polyhedron `P`.

# Example # unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
julia> P = V,FV;
julia> INERTIAMOMENT(P)
3-element PointNd:
 0.033333333333333326
 0.03333333333333334
 0.03333333333333334
```
# Example 2 # standard unit cube
```julia
V = [0.0 1.0 0.0 0.0 1.0 1.0 0.0 1.0; 0.0 0.0 1.0 0.0 1.0 0.0 1.0 1.0; 0.0 0.0 0.0 1.0 0.0 1.0 1.0 1.0]
FV= [[1,2,6],[6,4,1],[1,3,5],[5,2,1],[4,7,3],[3,1,4],[2,5,8],
[8,6,2],[3,7,8],[8,5,3],[4,6,8],[8,7,4]]
cube = V,FV
INERTIAMOMENT(cube)
3-element PointNd:
 0.6666666666666665
 0.6666666666666665
 0.6666666666666665
```
"""
function INERTIAMOMENT(P)::Array{Float64,1}
    out = zeros(3)
    result = SECONDMOMENT(P)
    out[1] = result[2] + result[3]
    out[2] = result[3] + result[1]
    out[3] = result[1] + result[2]
    return out
end



#///////////////////////////////////////////////////////////////
function chainAreas(V::Array{Float64,2},EV::Array{Int64,2},chains::Array{Int64,2})
	FE = [chains[:,f] for f=1:size(chains,2)]
	return chainAreas(V,EV,FE)
end

#///////////////////////////////////////////////////////////////
""" Implementation using integr.jl """
function chainAreas(V::Array{Float64,2}, EV::Array{Int64,2}, 
				chains::Array{Array{Int64,1},1})
	if size(V,1) == 2
		V = vcat(V,zeros(1,size(V,2)))
	end	
	pivots = [EV[:,abs(chain[1])][1] for chain in chains]
	out = zeros(length(pivots))
	for k=1:length(chains)
		area = 0
		triangles = [[] for h=1:length(chains[k])]
		for h=1:length(chains[k])
			edge = chains[k][h]
			v1,v2 = EV[:,abs(edge)]
			if sign(edge) == -1
				v1,v2=v2,v1
			end
			triangles[h] = Int[pivots[k],v1,v2]
		end
		P = V,hcat(triangles...)
		out[k] = Surface(P,true)
	end
	return out
end


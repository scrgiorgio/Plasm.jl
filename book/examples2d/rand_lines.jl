using Plasm

function randlines(n=400, t=0.4)
	V = zeros(Float64,2,2*n)
	EV = [zeros(Int64,2) for k=1:n]
	for k=1:n
		v1 = rand(Float64,2)
		v2 = rand(Float64,2)
		vm = (v1+v2)/2
		transl = rand(Float64,2)
		V[:,k] = (v1-vm)*t + transl
		V[:,n+k] = (v2-vm)*t + transl
		EV[k] = [k,n+k]
	end
	V,FVs,EVs = Lar.arrange2D(V,EV)
end

# ////////////////////////////////////////////////////////////
V,FVs,EVs = randlines();

# native OpenGL visualization
Plasm.VIEW(EXPLODER(V,FVs))

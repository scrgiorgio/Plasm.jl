

# //////////////////////////////////////////////////////////////////////
function constrained_triangulation2D(V::Points, EV::Cells)
	triin = Triangulate.TriangulateIO()
	triin.pointlist = V # scrgiorgio: this is by-col representation as LAR
	triin.segmentlist = hcat(EV...)
	(triout, __vorout) = Triangulate.triangulate("pQ", triin)  # exec triangulation
	return Array{Int64,1}[c[:] for c in eachcol(triout.trianglelist)]
end

# //////////////////////////////////////////////////////////////////////
function TRIANGULATE2D(V::Points, EV::Cells)
	num_vertices=size(V,2)
	copEV = lar2cop(EV)
	V_row = BYROW(V)
	trias = constrained_triangulation2D(V, EV)
	ret = Array{Int,1}[]
	for (u, v, w) in trias
		centroid = (V_row[u, :] + V_row[v, :] + V_row[w, :]) ./ 3
		if point_in_face(centroid, V_row, copEV)
			push!(ret, [u, v, w])
		end
	end
	return ret
end
export TRIANGULATE2D


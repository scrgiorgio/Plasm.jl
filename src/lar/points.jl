# LAR uses always by-col (!!!) representation
# if you need by-row use `V_row=BYROW(V)``
const Points = Matrix{Float64}
export Points


""" returns point dim """
function PDIM(V::Points)
	return size(V,1)
end
export PDIM

""" returns number of Lar vertices """
function NVERS(V::Points)
	return size(V,2)
end
export NVERS

# convert a vertex matrix from by-col (LAR default) to by-col
# assumption: input is by-col
function BYROW(V::Points)::Points
	return permutedims(V)
end
export BYROW

# convert a vertex matrix from by-row to by-col (LAR default)
# assumption: input is by-row
function BYCOL(V::Points)::Points
	return permutedims(V)
end
export BYCOL


# //////////////////////////////////////////////////////////////////////////////
""" Predicate to check equality of two vertices (only used above) """
function vertex_fuzzy_equals(v1, v2;err = 10e-8)
	return length(v1) == length(v2) && all(map((x1, x2) -> -err < x1 - x2 < err, v1, v2))
end
export vertex_fuzzy_equals

""" Predicate to check membership of `vertex` in `vertices_set` array"""
function is_visited_vertex(vertex, vertices_set)::Bool
	for v in vertices_set
		if vertex_fuzzy_equals(vertex, v)
			return true
		end
	end
	return false
end
export is_visited_vertex
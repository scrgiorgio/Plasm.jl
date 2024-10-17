# LAR uses always by-col (!!!) representation
# if you need by-row use `V_row=BYROW(V)` but try to avoid mixing different representations

const Points = Matrix{Float64}
export Points

function ToPoints(v::AbstractPointsNd)::Matrix{Float64}
	return hcat(to_concrete(v)...) 
end
export ToPoints


# /////////////////////////////////////////////////////////////
PointsDB=Dict{PointNd,Int}
export PointsDB

# ///////////////////////////////////////////////////////////////////
function round_vector(v::PointNd; digits::Int)::PointNd
  if digits==0 return v end
  return [round(value, digits=digits) for value in v]
end
export round_vector

# ///////////////////////////////////////////////////////////////////
function add_point(db::PointsDB, p::PointNd)::Int

  p=[abs(it)==0.0 ? 0.0 : it for it in p]

  if !haskey(db, p)  
    db[p]=length(db)+1 
  end
  return db[p]
end
export add_point

# ///////////////////////////////////////////////////////////////////
function get_points(db::PointsDB)::Points
	v=[PointNd() for I in 1:length(db)]
	for (pos,idx) in db
		v[idx]=pos
	end
	return hcat(v...)
end
export get_points


# /////////////////////////////////////////////////////////////////////
function point_distance(a, b)::Float64
  return LinearAlgebra.norm(a-b)
end

# /////////////////////////////////////////////////////////////////////
function point_is_between(a, b, c, epsilon):Bool
  return abs(point_distance(a,c) + point_distance(c,b) - point_distance(a,b)) <= epsilon
end

# /////////////////////////////////////////////////////////////////////
#  see https://gist.github.com/kylemcdonald/6132fc1c29fd3767691442ba4bc84018
function segment_intersection(p1, p2, p3, p4, epsilon::Float64)

  (x1,y1),(x2,y2),(x3,y3),(x4,y4) = p1,p2,p3,p4
  
  # parallel
  denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1)
  if abs(denom) < epsilon
    return nothing
  end

  # out of range
  ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom
  if (ua < -epsilon) || (ua > 1+epsilon) return nothing end

  # out of range
  ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom
  if (ub < -epsilon) || (ub > 1+epsilon) return nothing end

  return (
    x1 + ua * (x2-x1),
    y1 + ua * (y2-y1))
end

# ///////////////////////////////////////////////////////////////////
function bbox_create(points::Points)
  return (
    [minimum(points[R,:]) for R in 1:size(points,1)],
    [maximum(points[R,:]) for R in 1:size(points,1)]
  )
end
export bbox_create

# ///////////////////////////////////////////////////////////////////
function bbox_create(v::AbstractPointsNd)
  return bbox_create(ToPoints(v))
end

# ///////////////////////////////////////////////////////////////////
function bbox_intersect(box_a, box_b)
  (A,B),(C,D)=box_a,box_b
  return all([a <= d && b >= c for (a,b,c,d) in zip(A,B,C,D)])
end
export bbox_intersect

# ///////////////////////////////////////////////////////////////////
function bbox_contain(box_a, box_b)
  (A,B),(C,D)=box_a,box_b
  return all([a <= c && b >= d for (a,b,c,d) in zip(A,B,C,D)])
end
export bbox_contain

""" returns point dim (assuming by-col rep)"""
function PDIM(V::Points)
	return size(V,1)
end
export PDIM

""" returns number of Lar vertices (assuming by-col rep) """
function NVERTS(V::Points)
	return size(V,2)
end
export NVERTS

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

""" Predicate to check equality of two vertices (only used above ???) """
function vertex_fuzzy_equals(v1, v2; err = LAR_DEFAULT_ERR)
	return length(v1) == length(v2) && all([-err<(x-y)<err for (x,y) in zip(v1,v2)])
end
export vertex_fuzzy_equals


""" Predicate to check membership of `vertex` in `vertices_set` array"""
# posso solo generare vertici random
function is_visited_vertex(vertex, vertices_set)::Bool
	for v in vertices_set
		if vertex_fuzzy_equals(vertex, v)
			return true
		end
	end
	return false
end
export is_visited_vertex

# ////////////////////////////////////////////////////////////////////////
""" assuming by col """
function compute_centroid(V::Points)
  return Statistics.mean(eachcol(V))
end
export compute_centroid

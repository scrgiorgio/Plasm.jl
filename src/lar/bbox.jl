using IntervalTrees


"""
boundingbox(vertices::Points)

The axis aligned *bounding box* of the provided Matrix of n-dim `vertices`.
The box is returned as the pair of `Points` of two opposite corners.

NOTE: assuming LAR by-col representation, if it's by-row using dims=1
"""
function boundingbox(vertices::Points;dims::Int=2)
	minimum = mapslices(x -> min(x...), vertices, dims=dims)
	maximum = mapslices(x -> max(x...), vertices, dims=dims)
	return minimum, maximum
end
export boundingbox

function bbox_contains(container, contained)
	b1_min, b1_max = container
	b2_min, b2_max = contained
	all(map((i, j, k, l) -> i <= j <= k <= l, b1_min, b2_min, b2_max, b1_max))
end
export bbox_contains

function bbox_covering(bboxes, index, tree)
	covers = [[] for k = 1:length(bboxes)]
	for (i, boundingbox) in enumerate(bboxes)
		extent = bboxes[i][index, :]
		iterator = IntervalTrees.intersect(tree, tuple(extent...))
		for x in iterator
			append!(covers[i], x.value)
		end
	end
	return covers
end
export bbox_covering

""" Make dictionary of 1D boxes for IntervalTrees construction """
function bbox_coord_intervals(coord, bboxes)
	boxdict = OrderedDict{Array{Float64,1},Array{Int64,1}}()
	for (h, box) in enumerate(bboxes)
		key = box[coord, :]
		if haskey(boxdict, key) == false
			boxdict[key] = [h]
		else
			push!(boxdict[key], h)
		end
	end
	return boxdict
end
export bbox_coord_intervals

# //////////////////////////////////////////////////////////////////////////////
function bbox_containment_graph(bboxes)
	n = length(bboxes)
	ret = spzeros(Int8, n, n)
	for i in 1:n
		for j in 1:n
			if i != j && bbox_contains(bboxes[j], bboxes[i])
				ret[i, j] = 1
			end
		end
	end
	return ret
end
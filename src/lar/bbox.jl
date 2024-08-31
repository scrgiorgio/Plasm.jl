using IntervalTrees


# //////////////////////////////////////////////////////////////////////////////
"""
bbox(vertices::Points)

The axis aligned *bounding box* of the provided Matrix of n-dim `vertices`.
The box is returned as the pair of `Points` of two opposite corners.
"""
function bbox(vertices::Points)
	minimum = mapslices(x -> min(x...), vertices, dims=1)
	maximum = mapslices(x -> max(x...), vertices, dims=1)
	minimum, maximum
end
export bbox

function boundingbox(vertices::Points)
	minimum = mapslices(x -> min(x...), vertices, dims=2)
	maximum = mapslices(x -> max(x...), vertices, dims=2)
	return minimum, maximum
end
export boundingbox

function bbox_contains(container, contained)
	b1_min, b1_max = container
	b2_min, b2_max = contained
	all(map((i, j, k, l) -> i <= j <= k <= l, b1_min, b2_min, b2_max, b1_max))
end
export bbox_contains

function boxcovering(bboxes, index, tree)
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
export boxcovering

""" Make dictionary of 1D boxes for IntervalTrees construction """
function coordintervals(coord, bboxes)
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
export coordintervals
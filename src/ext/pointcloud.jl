# scrgiorgio 14Nov26: inserted for bpa and pointclouds
import PointClouds

# ///////////////////////////////////////////////////////////////////////
# read a LAS/LAZ file, returns a 3xN matrix (one point per column)
# if the file has more than `max_points` points, it is uniformly subsampled
function read_laz(filename::String; max_points::Int=typemax(Int))::Points
	las  = PointClouds.LAS(filename)
	N    = length(las)
	step = max(1, cld(N, max_points))
	# raw coordinates, no CRS transform (PointCloud() parses the CRS and fails on some files)
	coords = PointClouds.coordinates(step==1 ? las : las[1:step:N])
	ret = Points(undef, 3, length(coords))
	for (I, (x, y, z)) in enumerate(coords)
		ret[1, I], ret[2, I], ret[3, I] = x, y, z
	end
	return ret
end
export read_laz

# ///////////////////////////////////////////////////////////////////////
function cloud_geometry(V::Points, colors)::Geometry
	g = Geometry()
	g.cloud = copy(V)
	if !isnothing(colors)
		g.cloud_colors = colors
	end
	return g
end

# ///////////////////////////////////////////////////////////////////////
# build an Hpc storing V (one point per column) as a point cloud: no hulls, no topology
# `colors` is optional: RGBA per point, either a flat vector (4*N) or a 4xN matrix
# NOTE: unlike MKPOINTS (one hull per point) this scales to millions of points
#
# octree (when `octree_depth` > 0 or `octree_box` is given):
#   - `octree_box` = root 3D box: (min, max) corners e.g. ([0,0,0], [1,1,1]), or a Box3d / BoxNd;
#     default the point cloud bounding box
#   - points outside `octree_box` are dropped (with a warning)
#   - every node has 8 children, X varies fastest: child = ix + 2*iy + 4*iz + 1 (ix,iy,iz in 0:1)
#   - empty cells are kept (empty leaf), so the tree is always complete (8^depth leaves)
#   - leaves at `octree_depth` hold the points; cells are half-open [min,max) except the last one
#     along each axis that includes the box max, so every point goes into exactly one leaf
#   - each node has properties "octree_box" => (min, max) and "octree_depth" => level (root = 0)
#   - every node is an Hpc; a leaf has one child, the Geometry with the points (leaf.childs[1].cloud)
function POINTCLOUD(V::Points; colors=nothing, point_size::Int=DEFAULT_POINT_SIZE, point_color=DEFAULT_POINT_COLOR, octree_depth::Int=0, octree_box=nothing)::Hpc
	if !isnothing(colors)
		colors = Vector{Float32}(vec(colors))
		@assert length(colors) == 4*size(V, 2) "colors must have 4 values (RGBA) per point"
	end
	properties = Properties("point_size" => point_size, "point_color" => point_color)

	# no octree: a single geometry with all the points
	if octree_depth == 0 && isnothing(octree_box)
		return Hpc(MatrixNd(), [cloud_geometry(V, colors)], properties)
	end

	@assert octree_depth >= 0 "octree_depth must be >= 0"
	@assert size(V, 1) == 3 "the octree needs 3D points"

	# root box
	if isnothing(octree_box)
		p1, p2 = vec(minimum(V, dims=2)), vec(maximum(V, dims=2))
	else
		corners = octree_box isa Union{Box3d, BoxNd} ? (octree_box.p1, octree_box.p2) : octree_box
		p1, p2 = Vector{Float64}(collect(corners[1])), Vector{Float64}(collect(corners[2]))
		@assert length(p1) == 3 && length(p2) == 3 && all(p1 .<= p2) "octree_box must be (min, max) with min <= max"
		inside = vec(all(p1 .<= V .<= p2, dims=1))
		dropped = count(!, inside)
		if dropped > 0
			@warn "POINTCLOUD: $dropped of $(size(V, 2)) points are outside octree_box $((p1, p2)) and are dropped"
			V = V[:, inside]
			colors = isnothing(colors) ? nothing : colors[repeat(inside, inner=4)]
		end
	end

	# leaf cell of each point (0-based), grouped by leaf
	n = 2^octree_depth
	box_size = p2 .- p1
	leaf_key(c) = c[1] + n * c[2] + n * n * c[3]
	keys = Vector{Int}(undef, size(V, 2))
	for I in 1:size(V, 2)
		c = [box_size[A] > 0 ? clamp(floor(Int, (V[A, I] - p1[A]) / box_size[A] * n), 0, n - 1) : 0 for A in 1:3]
		keys[I] = leaf_key(c)
	end
	order = sortperm(keys)
	leaf_points = Dict{Int,Vector{Int}}()
	for I in order
		push!(get!(leaf_points, keys[I], Int[]), I)
	end

	# build the tree; c = cell coordinates (0-based) at `level`
	function build(level::Int, c::Vector{Int})::Hpc
		cell_size = box_size ./ 2^level
		node_properties = Properties("octree_box" => (p1 .+ c .* cell_size, p1 .+ (c .+ 1) .* cell_size), "octree_depth" => level)
		if level == octree_depth
			idx = get(leaf_points, leaf_key(c), Int[])
			leaf_colors = isnothing(colors) ? nothing : vec(reshape(colors, 4, :)[:, idx])
			return Hpc(MatrixNd(), [cloud_geometry(V[:, idx], leaf_colors)], node_properties)
		end
		childs = [build(level + 1, 2 .* c .+ [ix, iy, iz]) for iz in 0:1 for iy in 0:1 for ix in 0:1] # X fastest
		return Hpc(MatrixNd(), childs, node_properties)
	end

	root = build(0, [0, 0, 0])
	merge!(root.properties, properties)
	return root
end
export POINTCLOUD

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
# build an Hpc storing V (one point per column) as a point cloud: no hulls, no topology
# `colors` is optional: RGBA per point, either a flat vector (4*N) or a 4xN matrix
# NOTE: unlike MKPOINTS (one hull per point) this scales to millions of points
function POINTCLOUD(V::Points; colors=nothing, point_size::Int=DEFAULT_POINT_SIZE, point_color=DEFAULT_POINT_COLOR)::Hpc
	g = Geometry()
	g.cloud = copy(V)
	if !isnothing(colors)
		g.cloud_colors = Vector{Float32}(vec(colors))
		@assert length(g.cloud_colors) == 4*size(V, 2) "colors must have 4 values (RGBA) per point"
	end
	return Hpc(MatrixNd(), [g], Properties("point_size" => point_size, "point_color" => point_color))
end
export POINTCLOUD

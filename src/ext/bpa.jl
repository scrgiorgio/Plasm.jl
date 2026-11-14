# scrgiorgio 14Nov26: inserted for bpa and pointclouds
import BPA as BallPivoting

# ///////////////////////////////////////////////////////////////////////
# all point cloud points of an Hpc (transformations applied), 3xN by-col
function cloud_points(obj::Hpc)::Points
	ret = Points(undef, 3, 0)
	for (T, properties, g) in toList(obj)
		isempty(g.cloud) && continue
		D, N = size(g.cloud)
		Tdim = max(dim(T), D + 1)
		H = embed(T, Tdim).T * vcat(ones(1, N), g.cloud, zeros(Tdim - D - 1, N)) # homogeneous coordinate is the first
		P = H[2:end, :] ./ H[1:1, :]
		ret = hcat(ret, vcat(P[1:min(3, end), :], zeros(max(0, 3 - size(P, 1)), N)))
	end
	return ret
end
export cloud_points

# ///////////////////////////////////////////////////////////////////////
# typical distance between neighbouring points (median nearest-neighbour distance, on a random sample)
function POINTSPACING(obj::Hpc)::Float64
	return BallPivoting.estimate_spacing([BallPivoting.Vec3(p...) for p in eachcol(cloud_points(obj))])
end
export POINTSPACING

# ///////////////////////////////////////////////////////////////////////
# Ball Pivoting surface reconstruction (https://github.com/ctsilva/BPA.jl)
#   obj      Hpc with point clouds (see POINTCLOUD)
#   radii    ball radius (or radii, one pass each); default 2x and 4x the point spacing
#   normals  3xN outward normals; default estimated from the `k` nearest neighbours
#   up       if given (e.g. [0,0,1] for aerial LIDAR) estimated normals are flipped to point towards it
# returns an Hpc with one triangle per hull
function BALLPIVOTING(obj::Hpc; radii=nothing, normals=nothing, k::Int=10, up=nothing)::Hpc
	V = cloud_points(obj)
	@assert size(V, 2) >= 3 "BALLPIVOTING needs an Hpc with a point cloud"
	positions = [BallPivoting.Vec3(p...) for p in eachcol(V)]

	if isnothing(normals)
		normals = BallPivoting.estimate_normals(positions; k=k)
		if !isnothing(up)
			normals = [n' * up < 0 ? -n : n for n in normals]
		end
	else
		normals = [BallPivoting.Vec3(n...) for n in eachcol(normals)]
	end

	if isnothing(radii)
		spacing = BallPivoting.estimate_spacing(positions)
		radii = [2.0, 4.0] .* spacing
	end

	mesh = BallPivoting.reconstruct(BallPivoting.PointCloud(positions, normals), radii isa Real ? [radii] : radii)
	points = [V[:, I] for I in 1:size(V, 2)]
	return Hpc(CreateGeometry(points, [collect(t) for t in mesh.triangles]))
end
export BALLPIVOTING

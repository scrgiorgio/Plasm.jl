# Read a LAZ file and render it as a point cloud, optionally stored in an octree (usage and options: see README.md)

using Plasm

# command line: one file name plus --name=value options
function option(name, default)
	values = [arg[length(name)+4:end] for arg in ARGS if startswith(arg, "--$name=")]
	return isempty(values) ? default : last(values)
end

files = filter(arg -> !startswith(arg, "--"), ARGS)
length(files) == 1 || error("usage: julia --project=. tutorials/point-cloud.jl <file.laz> [--color-by-height] [--max-points=2000000] [--octree-depth=0] [--octree-box=xmin,xmax,ymin,ymax,zmin,zmax]")
filename        = files[1]

color_by_height = "--color-by-height" in ARGS

max_points      = parse(Int, option("max-points", "2000000"))

octree_depth    = parse(Int, option("octree-depth", "0"))

# octree box in file coordinates, min and max for each axis (default: point cloud bounding box)
octree_box      = option("octree-box", nothing)
if !isnothing(octree_box)
	octree_box = parse.(Float64, Base.split(octree_box, ","))
	length(octree_box) == 6 || error("--octree-box needs 6 values: xmin,xmax,ymin,ymax,zmin,zmax")
end

V = read_laz(filename, max_points=max_points)
println("Read ", size(V, 2), " points")
p1, p2 = minimum(V, dims=2), maximum(V, dims=2)
println("Bounding box: --octree-box=", join(round.([p1[1], p2[1], p1[2], p2[2], p1[3], p2[3]], digits=3), ","))

# move to the origin (UTM coordinates are too big for the Float32 used by the viewer)
V = V .- p1
if !isnothing(octree_box)
	xmin, xmax, ymin, ymax, zmin, zmax = octree_box
	octree_box = ([xmin, ymin, zmin] .- p1[:], [xmax, ymax, zmax] .- p1[:])
end

# (optional) color points by elevation (blue=low, red=high)
colors = nothing
if color_by_height
	zmin, zmax = extrema(V[3, :])
	colors = Vector{Float32}()
	for z in V[3, :]
		t = (z - zmin) / (zmax - zmin)
		append!(colors, Point4d(t, 1.0 - abs(2t - 1.0), 1.0 - t, 1.0))
	end
end

# store the points in an Hpc (without colors, `point_color` is used for all points)
# with octree_depth > 0 the points go into an octree with 8^octree_depth leaves
@time obj = POINTCLOUD(V, 
	colors=colors, 
	point_size=octree_depth > 0 ? 2 : 1, 
	point_color=BLACK, 
	octree_depth=octree_depth,
	octree_box=octree_box)

# octree: give every leaf (the node holding the Geometry) its own color
function color_leaves(node)
	if isa(node.childs[1], Geometry)
		node.properties["point_color"] = Point4d(rand(), rand(), rand(), 1.0)
	else
		foreach(color_leaves, node.childs)
	end
end
if octree_depth > 0
	color_leaves(obj)
end

VIEW(obj,
	title=basename(filename),
	properties=Properties(
		"background_color" => Point4d(1.0, 1.0, 1.0, 1.0),
		"show_axis" => true
		)
)

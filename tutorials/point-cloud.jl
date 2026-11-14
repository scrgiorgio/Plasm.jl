# Read a LAZ file and render it as a point cloud (usage and options: see README.md)

using Plasm

# command line: one file name plus --name=value options
function option(name, default)
	values = [arg[length(name)+4:end] for arg in ARGS if startswith(arg, "--$name=")]
	return isempty(values) ? default : last(values)
end

files = filter(arg -> !startswith(arg, "--"), ARGS)
length(files) == 1 || error("usage: julia --project=. tutorials/point-cloud.jl <file.laz> [--color-by-height] [--max-points=2000000]")
filename        = files[1]

color_by_height = "--color-by-height" in ARGS

max_points      = parse(Int, option("max-points", "2000000"))

V = read_laz(filename, max_points=max_points)
println("Read ", size(V, 2), " points")

# move to the origin (UTM coordinates are too big for the Float32 used by the viewer)
V = V .- minimum(V, dims=2)

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
@time obj = POINTCLOUD(V, colors=colors, point_size=1, point_color=BLACK)

VIEW(obj, 
	title=basename(filename),
	properties=Properties(
		"background_color" => Point4d(1.0, 1.0, 1.0, 1.0),
		"show_axis" => true
		)
)

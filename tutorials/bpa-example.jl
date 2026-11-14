# Surface reconstruction of a LIDAR point cloud with the Ball Pivoting Algorithm (usage and options: see README.md)

using Plasm

# command line: one file name plus --name=value options
function option(name, default)
	values = [arg[length(name)+4:end] for arg in ARGS if startswith(arg, "--$name=")]
	return isempty(values) ? default : last(values)
end
files = filter(arg -> !startswith(arg, "--"), ARGS)
length(files) == 1 || error("usage: julia --project=. tutorials/bpa-example.jl <file.laz> [--radii=2,4] [--max-points=200000]")
filename   = files[1]
factors    = parse.(Float64, Base.split(option("radii", "2,4"), ","))
max_points = parse(Int, option("max-points", "200000"))

V = read_laz(filename, max_points=max_points)
println("Read ", size(V, 2), " points")

# move to the origin (UTM coordinates are too big for the Float32 used by the viewer)
V = V .- minimum(V, dims=2)

# point cloud stored in an Hpc
cloud = POINTCLOUD(V, point_size=1, point_color=BLACK)

# ball radii from the point spacing
spacing = POINTSPACING(cloud)
radii = factors .* spacing
println("Point spacing: ", round(spacing, digits=3), ", ball radii: ", round.(radii, digits=3), " (", factors, " x spacing)")

# reconstruct the surface: normals are estimated from the points, and since this is
# an aerial scan they are oriented upwards
@time mesh = BALLPIVOTING(cloud, radii=radii, up=[0.0, 0.0, 1.0])
used = length(mesh.childs[1].points)
percent = round(Int, 100 * used / size(V, 2))
println("Triangles: ", length(mesh.childs[1].hulls), ", points used: ", used, " of ", size(V, 2), " (", percent, "%)")

# tips for the next run: --radii are multiples of the spacing (radius = factor x spacing)
radii_str(f) = join([isinteger(x) ? string(Int(x)) : string(x) for x in round.(f, digits=2)], ",")
println("Tips (point spacing ", round(spacing, digits=3), "): a reasonable start is --radii=2,4 (radius ",
	round(2 * spacing, digits=2), " and ", round(4 * spacing, digits=2), ")")
if minimum(factors) < 1.5
	println("  - smallest radius is too close to the spacing: the ball falls through the points, try --radii=2,4")
elseif percent < 97 && maximum(factors) < 8
	println("  - ", percent, "% of the points used: to fill the holes add a larger pass, try --radii=", radii_str([factors; 2 * maximum(factors)]))
elseif percent < 97
	println("  - ", percent, "% of the points used: larger radii no longer help, the unused points are probably isolated (noise, sparse vegetation)")
else
	println("  - almost all points used (", percent, "%): if details look bridged over (trees, walls), try smaller radii --radii=", radii_str(max.(1.5, factors ./ 2)))
end

# input points on the left, reconstructed surface on the right (press L in the viewer to show the triangle edges)
width = maximum(V[1, :]) * 1.1
VIEW(
	STRUCT(
		cloud,
		T(1)(width)(PROPERTIES(mesh, Properties(
			"face_color" => Point4d(0.8, 0.8, 0.8, 1.0),
			"line_color" => Point4d(0.3, 0.3, 0.3, 1.0)
		)))
	),
	title=basename(filename),
	properties=Properties(
		"background_color" => Point4d(1.0, 1.0, 1.0, 1.0),
		"show_axis" => false,
		"show_lines" => false
	)
)

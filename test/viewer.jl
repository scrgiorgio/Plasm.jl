using Plasm

function MyMain()

	GLView([
		GLCuboid(Box3d(Point3d(0,0,0),Point3d(1,1,1)))
		GLAxis(Point3d(0,0,0),Point3d(+1.1,+1.1,+1.1))
		])	

	GLView(GLText("hello"))


	# example of customizing the viewer position
	if false
		VIEW(
			CUBOID([1,1,1]),
			Dict(
				# 3d pipeline position
				"pos" => Point3d(0.5, 0.5, 3.0),
				"dir" => Point3d(0.0, 0.0,-1.0),
				"vup" => Point3d(1.0, 0.0, 0.0),
				"znear" => 0.1,
				"zfar"  => 10.0,
				# perspective fov
				"fov" => 60.0,
				#triangles, show/hide lines
				"show_lines" => false,
				#viewer background color
				"background_color" => Point4d(1.0, 1.0, 1.0, 1.0),
				# perspective or ortho projection
				"use_ortho" => true
			)
		)
	end
	
	
	begin
		
		V = [9. 13 15 17 14 13 11 9 7 5 3 0 2 2 5 7 4 12 6 8 3 5 7 8 10 11 10 13 14 13 11 9 7 4 2 12 12; 0 2 4 8 9 10 11 10 9 9 8 6 3 1 0 1 2 10 3 3 5 5 6 5 5 4 2 4 6 7 9 7 7 7 6 7 5];
		EV = Array{Int64,1}[[1, 2], [1, 16], [1, 27], [2, 3], [2, 27], [2, 28], [3, 4], [3, 29], [4, 5], [4, 29], [5, 6], [5, 29], [5, 30], [6, 7], [6, 18], [7, 8], [7, 18], [8, 9], [8, 31], [8, 32], [9, 10], [9, 33], [10, 11], [10, 34], [11, 12], [11, 34], [12, 13], [12, 21], [12, 35], [13, 14], [13, 17], [13, 21], [14, 15], [14, 17], [15, 16], [15, 17], [16, 20], [16, 27], [17, 19], [18, 31], [19, 20], [19, 22], [19, 23], [19, 24], [20, 24], [20, 25], [21, 22], [21, 34], [21, 35], [22, 23], [22, 33], [23, 24], [23, 33], [24, 25], [24, 32], [25, 26], [25, 27], [25, 32], [25, 36], [26, 27], [26, 28], [26, 37], [28, 29], [28, 37], [29, 30], [29, 37], [30, 36], [31, 32], [31, 36], [32, 33], [33, 34], [34, 35], [36, 37]];
	
		obj=Hpc(V,EV)
	
		batches=Vector{GLBatch}()
	
		append!(batches,GetBatchesForHpc(obj))
	
		append!(batches,GLText(
			[V[:,k] for k=1:size(V,2)],
			EV=[it for it in EV],
			# FV=FV,
			V_color=Point4d(1,1,1,1),
			EV_color=Point4d(1,0,0,1),
			FV_color=Point4d(0,1,0,1)))
	
		View(batches)
	
	end

	
end

MyMain()
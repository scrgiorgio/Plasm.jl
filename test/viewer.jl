using Plasm

# ////////////////////////////////////////////////////////////////////////
function GLTriangle(p1,p2,p3, n)

	triangles=GLBatch(TRIANGLES)

	append!(triangles.vertices.vector,p1)
	append!(triangles.vertices.vector,p2)
	append!(triangles.vertices.vector,p3)

	append!(triangles.normals.vector,n)
	append!(triangles.normals.vector,n)
	append!(triangles.normals.vector,n)

	return triangles
end

function MyMain()

	if false
		n_up=Point3d(0.0, 0.0, +1.0)
		n_dw=Point3d(0.0, 0.0, -1.0)

		GLView([
			GLTriangle(Point3d(0.0, 0.0, 0.0),Point3d(1.0, 0.0, 0.0),Point3d(0.0 ,1.0, 0.0), n_up)
			GLTriangle(Point3d(0.0, 0.0, 1.0),Point3d(1.0, 0.0, 1.0),Point3d(0.0 ,1.0, 1.0), n_dw)
			GLTriangle(Point3d(0.0, 0.0, 2.0),Point3d(0.0, 1.0, 2.0),Point3d(1.0 ,0.0, 2.0), n_up)
			GLTriangle(Point3d(0.0, 0.0, 3.0),Point3d(0.0, 1.0, 3.0),Point3d(1.0, 0.0, 3.0), n_dw)

			GLAxis(Point3d(0,0,0),Point3d(+1.1,+1.1,+1.1))
			])	
	else
	
		GLView([
			GLCuboid(Box3d(Point3d(0,0,0),Point3d(1,1,1)))
			GLAxis(Point3d(0,0,0),Point3d(+1.1,+1.1,+1.1))
			])	
	end
	
end

MyMain()
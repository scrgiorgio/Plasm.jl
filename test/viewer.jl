using Plasm

function MyMain()

	GLView([
		GLCuboid(Box3d(Point3d(0,0,0),Point3d(1,1,1)))
		GLAxis(Point3d(0,0,0),Point3d(+1.1,+1.1,+1.1))
		])	

	GLView([GLText("hello")])
	
end

MyMain()
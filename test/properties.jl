using Plasm
VIEW(CUBE(1),title="normal cube (1)")

VIEW(CUBE(1),"normal cube (2)")
VIEW(CUBE(1),title="normal cube (13)")
VIEW(CUBE(1),background_color=[0.8,0.1,0.5],title="cube with background color")

VIEW(PROPERTIES(CUBE(1),Dict("line_width"=>10)),title="line_width")
VIEW(PROPERTIES(CUBE(1),Dict("line_color"=>[1.0,0.0,0.0])),title="line_color")
VIEW(PROPERTIES(CUBE(1),Dict("face_color"=>[0.0,1.0,0.0])),title="face_color")


VIEW(STRUCT(
	PROPERTIES(STRUCT(T(1)(0.0),CUBE(1)), Dict("face_color"=>[0.0,1.0,0.0])),
	PROPERTIES(STRUCT(T(1)(1.0),CUBE(1)), Dict("face_color"=>[1.0,0.0,0.0]))),
	title="2 colored cubes"
)


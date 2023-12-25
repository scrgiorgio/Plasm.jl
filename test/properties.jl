using Plasm

# different ways to `View`
VIEW(CUBE(1),title="normal cube (1)")
VIEW(CUBE(1), "normal cube (2)")
VIEW(CUBE(1),title="normal cube (13)")

# example: how to set the overall background color
VIEW(CUBE(1),background_color=[0.8,0.1,0.5],title="cube with background color")

# line_with property
VIEW(
   PROPERTIES(
      CUBE(1),
      Dict("line_width"=>10)),
   title="line_width")

# line_color property
VIEW(
  PROPERTIES(
     CUBE(1),
     Dict("line_color"=>Point4d(1.0,0.0,0.0,1.0))),
  title="line_color"
)

# face_color property
VIEW(
   PROPERTIES(
     CUBE(1),
     Dict("face_color"=>Point4d(0.0,1.0,0.0,1.0))),
   title="face_color"
)

# example of mixing and matching different properties
VIEW(STRUCT(
   PROPERTIES(
      STRUCT(T(1)(0.0),CUBE(1)), 
      Dict(
         "face_color"=>Point4d(0.0,1.0,0.0,1.0),
         "line_color"=>Point4d(1.0,1.0,0.0,1.0),
         "line_width"=>3
      )
   ),
   PROPERTIES(
      STRUCT(T(1)(1.0),CUBE(1)), 
      Dict(
         "face_color"=>Point4d(1.0,0.0,0.0,1.0),
         "line_color"=>Point4d(0.0,1.0,1.0,1.0),
         "line_width"=>3
      )
   )),
   title="2 colored cubes",
   background_color=[0.0,0.0,0.0]
)

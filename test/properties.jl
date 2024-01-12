using Plasm

# different ways to `View`
VIEW(CUBE(1), "normal cube (1)")

# example: how to set the overall background color
VIEW(CUBE(1),
   Dict(
      "background-color" => [0.8,0.1,0.5],
      "title" => "cube with background color"
   ))

# line_with property
VIEW(
   PROPERTIES(
      CUBE(1),
      Dict("line_width"=>10)),
   "line_width")

# line_color property
VIEW(
  PROPERTIES(
     CUBE(1),
     Dict("line_color"=>Point4d(1.0,0.0,0.0,1.0))),
  "line_color"
)

# face_color property
VIEW(
   PROPERTIES(
     CUBE(1),
     Dict("face_color"=>Point4d(0.0,1.0,0.0,1.0))),
   "face_color"
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

   Dict(
      "title" => "2 colored cubes",
      "background-color" => [0.0,0.0,0.0]
   )
)

# //////////////////////////////////////////////////////
# 2d frame/hull/point example
begin

   points=[[rand()+1.0,rand()+1.0] for I in 1:100]
 
   obj=STRUCT(
 
     # show points
     PROPERTIES(
       MKPOINTS(points),
       Dict(
         "point_color"=>YELLOW, 
         "point_size"=>3
       )
     ),
     
     # show hull (note face color is transparent, so no fill)
     PROPERTIES(
       MKPOL(
         points,
         [[it for it in 1:length(points)]]
       ),
       Dict(
         "face_color"=>TRANSPARENT,
         "line_color"=>GREEN, 
         "line_width"=>2
       )
     ),
     # show frame
     FRAME2([0.0,0.0],[1.0,1.0]),
   )
 
   VIEW(obj, Dict("show-axis" => false))
end
 
# //////////////////////////////////////////////////////
# 3d frame/hull/point example
begin

   points=[[rand()+1.0,rand()+1.0,rand()+1.0] for I in 1:100]

   obj=STRUCT(

      # show points
      PROPERTIES(
         MKPOINTS(points),
         Dict(
         "point_color"=>YELLOW, 
         "point_size"=>3
         )
      ),
      # show hull
      PROPERTIES(
         MKPOL(
         points,
         [[it for it in 1:length(points)]]
         ),
         Dict(
         "face_color"=>GRAY,
         "line_color"=>GREEN, 
         "line_width"=>2
         )
      ),
      # show frame
      FRAME3(Point3d(0,0,0),Point3d(1,1,1)),
   )

   VIEW(obj, Dict("show-axis" => false) )

end
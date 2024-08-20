using Plasm

"""
# List of properties
View(batches, properties=Properties(
  "background_color" => Point4d(1,1,1,1))
  "title"            => "my title")
  "use_ortho"        => false)
  "show_lines"       => true)
  "fov"              => 60.0)
  "pos"              => Point3d(0,0,0) )
  "dir"              => Point3d(0,0,-1))
  "vup"              => Point3d(1,0,0))
  "znear"            => 0.00001) 
  "zfar "            => 100.0) 
  "walk_speed "      => 10.0) 
)
"""

# different ways to `View`
VIEW(CUBE(1), "normal cube (1)")

# example: how to set the overall background color
VIEW(CUBE(1),
Properties(
      "background_color" => [0.8,0.1,0.5],
      "title" => "cube with background color"
   ))

# line_with property
VIEW(
   PROPERTIES(
      CUBE(1),
      Properties("line_width"=>10)),
   "line_width")

# line_color property
VIEW(
  PROPERTIES(
     CUBE(1),
     Properties("line_color"=>Point4d(1.0,0.0,0.0,1.0))),
  "line_color"
)

# face_color property
VIEW(
   PROPERTIES(
     CUBE(1),
     Properties("face_color"=>Point4d(0.0,1.0,0.0,1.0))),
   "face_color"
)

# COLOR
VIEW(COLOR([1.0,1.0,0.0,0.0])(CUBE(1)))

# example of mixing and matching different properties
cube_a=STRUCT(T(1)(0.0),CUBE(1))
a=PROPERTIES(cube_a, 
Properties(
      "face_color"=>Point4d(0.0,1.0,0.0,1.0),
      "line_color"=>Point4d(1.0,1.0,0.0,1.0),
      "line_width"=>3
   )
)
cube_b=STRUCT(T(1)(1.0),CUBE(1))
b=PROPERTIES(cube_b, 
Properties(
      "face_color"=>Point4d(1.0,0.0,0.0,1.0),
      "line_color"=>Point4d(0.0,1.0,1.0,1.0),
      "line_width"=>3
   )
)
VIEW(
   STRUCT(a,b),
   Properties(
      "title" => "2 colored cubes",
      "background_color" => [0.0,0.0,0.0]
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
       Properties(
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
       Properties(
         "face_color"=>TRANSPARENT,
         "line_color"=>GREEN, 
         "line_width"=>2
       )
     ),
     # show frame
     FRAME2([0.0,0.0],[1.0,1.0]),
   )
 
   VIEW(obj, Properties("show-axis" => false))
end
 
# //////////////////////////////////////////////////////
# 3d frame/hull/point example
begin

   points=[[rand()+1.0,rand()+1.0,rand()+1.0] for I in 1:100]

   obj=STRUCT(

      # show points
      PROPERTIES(
         MKPOINTS(points),
         Properties(
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
         Properties(
         "face_color"=>GRAY,
         "line_color"=>GREEN, 
         "line_width"=>2
         )
      ),
      # show frame
      FRAME3(Point3d(0,0,0),Point3d(1,1,1)),
   )

   VIEW(obj, Properties("show_axis" => false) )

end
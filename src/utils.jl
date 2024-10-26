
# /////////////////////////////////////////////////////////////////////
function print_matrix(name::String,value::Matrix)
  println(name)
  for (I,it) in enumerate(eachcol(value)) 
    println(I," ",it) 
  end
end

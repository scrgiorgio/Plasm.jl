
# /////////////////////////////////////////////////////////////////////
function print_matrix(name::String,value::Matrix)
  println(name)
  for (I,it) in enumerate(eachcol(value)) 
    println(I," ",it) 
  end
end

# /////////////////////////////////////////////////////////////////////
function get_number_of_digits(value::Float64) 
	@assert(value<=0.1)
	ret = 1 
	while round(value,digits=ret)!=value
		ret += 1 
	end
	return ret 
end

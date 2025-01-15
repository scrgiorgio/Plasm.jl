
using Plasm
using EzXML
using DataStructures

# ///////////////////////////////////////////////
function get_svg_class(node)
  try
    return node["class"]
catch
    return "default"
end
end

# ///////////////////////////////////////////////
function parse_svg(filename::String; flag=true)

  ret=Dict()
  body = read(filename,String)
  doc = parsexml(body)
  for it in eachelement(doc.root)

    # ___________________________________________
    if it.name=="line"
      class=get_svg_class(it)
      if !haskey(ret,class) ret[class]=[] end
      x1, y1, x2, y2=[parse(Float64, it[k]) for k in ("x1","y1","x2","y2")]
      println("line $(x1) $(y1) $(x2) $(y2)")
      push!(ret[class],[x1,y1, x2,y2])
      continue
    end

    # ___________________________________________
    if it.name=="rect"
      class=get_svg_class(it)
      if !haskey(ret,class) ret[class]=[] end
      x1, y1, width, height=[parse(Float64, it[k]) for k in ("x","y","width","height")]
      x2, y2 = x1+width, y1+height
      println("rect $(x1) $(y1) $(x2) $(y2)")
      push!(ret[class],[x1,y1, x2,y1])
      push!(ret[class],[x2,y1, x2,y2])
      push!(ret[class],[x2,y2, x1,y2])
      push!(ret[class],[x1,y2, x1,y1])
      continue
    end

    # ___________________________________________
    if it.name=="path"
      class=get_svg_class(it)
      if !haskey(ret,class) ret[class]=[] end
      # println("path")
      startx,starty=nothing,nothing # first point
      lastx,lasty=nothing,nothing

			d=replace(it["d"],"," => " ")

      tokens= [jt for jt in Base.split(d," ") if length(jt)>0] 
      
      I=1
      while I<=length(tokens)

        # eatc token
        cmd=tokens[I]
        I+=1

        # Close Path
        if cmd=="Z" || cmd=="z"
          # println("  Close path")
          push!(ret[class],[lastx,lasty, startx,starty])
          continue
        end

        numbers=[]
        while I<=length(tokens) && !(length(tokens[I])==1 && isletter(tokens[I][1]))
          push!(numbers, parse(Float64, tokens[I]) )
          I+=1
        end

        # Move To
        if cmd=="M" || cmd=="m"
          # println("  MoveTo ",numbers)
          lastx,lasty=numbers
          startx,starty=lastx,lasty

        # line to
        elseif cmd=="L" || cmd=="l"
          # println("  LineTo ",numbers)
          x2,y2=numbers
          push!(ret[class],[lastx,lasty, x2,y2])
          lastx,lasty=x2,y2

        # vertical line
				elseif cmd=="V" || cmd=="v"
					x2,y2=lastx, numbers[1]
					push!(ret[class],[lastx,lasty, x2,y2])
					lastx,lasty=x2,y2

        # horizontal line
				elseif cmd=="H" || cmd=="h"
					x2,y2=numbers[1], lasty
					push!(ret[class],[lastx,lasty, x2,y2])
					lastx,lasty=x2,y2

        else
          println("Unsupported cmd=[$(cmd)]")
          @assert(false) # unsupported
        end
      end

      continue

    end

    # not supported
    println("Skipping $(it.name)")
    continue
  end

  return ret

end



# ///////////////////////////////////////////////
function read_svg(filename::String; line_width::Int=1)::Hpc
  v=[]
	for (class,lines) in parse_svg(filename)

		println("/////////// class $(class)")

		points,hulls=Plasm.ConcretePointsNd(),Cells()
		for (x1,y1,x2,y2) in lines
      push!(points, [x1,y1])
      push!(points, [x2,y2])
      N=length(points)
      push!(hulls,[N-1,N])
		end
		push!(v,PROPERTIES(
			MKPOL(points, hulls),
			Properties(
        "line_color" => Plasm.COLORS[rand(1:12)], 
        "line_width" => line_width)
		))
	end

	return STRUCT(v)
end
export read_svg



"""

"""

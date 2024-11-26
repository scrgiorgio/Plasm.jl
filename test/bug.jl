using Plasm
using EzXML

# ///////////////////////////////////////////////
function split_any(s, chars)::Vector{String}
  sep=string(chars[1])
  for other in chars[2:end]
    s=string(replace(s, string(other) => sep))
  end
  return Base.split(s,sep)
end

# ///////////////////////////////////////////////
function parse_svg_path(d)::Vector{String}
  prev,ret=nothing,[]
  for I in eachindex(d)
    if isletter(d[I])
      if !isnothing(prev) 
        token=string(strip(d[prev:I-1]))
        if length(token)>0
          push!(ret,token) 
        end
      end
      push!(ret,string(d[I]))
      prev=I+1
    end
  end
  return ret
end

# //////////////////////////////////////////////////////////////////////////////
"""
	lines2lar(lines)

LAR model construction from array of float quadruples.
Each `line` in input array stands for `x1,y1,x2,y2`.
"""
function lines2lar(lines)
	vertdict = OrderedDict{Array{Float64,1}, Int64}()
	EV = Array{Int64,1}[]
	idx = 0
	for h=1:size(lines,2)
		x1,y1,x2,y2 = lines[:,h]

		if ! haskey(vertdict, [x1,y1])
			idx += 1
			vertdict[[x1,y1]] = idx
		end
		if ! haskey(vertdict, [x2,y2])
			idx += 1
			vertdict[[x2,y2]] = idx
		end
		v1,v2 = vertdict[[x1,y1]],vertdict[[x2,y2]]
		push!(EV, [v1,v2])
	end
	V = hcat(collect(keys(vertdict))...)
	return V,EV
end


# //////////////////////////////////////////////////////////////////////////////
"""
	svg_normalize(V::Points; flag=true::Bool)::Points

2D normalization transformation (isomorphic by defaults) of model
vertices to normalized coordinates ``[0,1]^2``. Used with SVG importing.
"""
function svg_normalize(V::Points; flag=true)
	m,n = size(V)
	if m > n # V by rows
		V = convert(Points, V')
	end

	xmin = minimum(V[1,:]); ymin = minimum(V[2,:]);
	xmax = maximum(V[1,:]); ymax = maximum(V[2,:]);
	box = [[xmin; ymin] [xmax; ymax]]	# containment box
	aspectratio = (xmax-xmin)/(ymax-ymin)
	if flag
		if aspectratio > 1
			umin = 0; umax = 1
			vmin = 0; vmax = 1/aspectratio; ty = vmax
		elseif aspectratio < 1
			umin = 0; umax = aspectratio
			vmin = 0; vmax = 1; ty = vmax
		end
		Z = T(1,2)(0,ty) ∘ S(1,2)(1,-1) ∘ S(1,2)(umax-umin,vmax-vmin) ∘
			S(1,2)(1/(xmax-xmin),1/(ymax-ymin)) ∘ T(1,2)(-xmin,-ymin)
	else
		# Z = T(1,2)(0, ymax-ymin) * S(1,2)(1,-1)
		Z = S(1,2)(1,-1)
	end
	return Z # Plasm tensor = composition of transformations
end


# ///////////////////////////////////////////////
function parse_svg(filename::String; flag=true)

  lines=[]
  body = read(filename,String)
  doc = parsexml(body)
  for it in eachelement(doc.root)

    # ___________________________________________
    if it.name=="line"
      x1, y1, x2, y2=[parse(Float64, it[k]) for k in ("x1","y1","x2","y2")]
      push!(lines,[x1,y1, x2,y2])
      continue
    end

    # ___________________________________________
    if it.name=="rect"
      x1, y1, width, height=[parse(Float64, it[k]) for k in ("x","y","width","height")]
      x2,y2=x1+width,y1+height
      push!(lines,[x1,y1, x2,y1])
      push!(lines,[x2,y1, x2,y2])
      push!(lines,[x2,y2, x1,y2])
      push!(lines,[x1,y2, x1,y1])
      continue
    end

    # ___________________________________________
    if it.name=="path"
      X1,Y1=nothing,nothing # first point
      x1,y1=nothing,nothing
      tokens=parse_svg_path(it["d"])
      I=1;while I<=length(tokens)

        cmd=tokens[I];I+=1

        # Move To
        if cmd=="M" || cmd=="m"
          x1,y1=[parse(Float64, jt) for jt in split_any(tokens[I]," ,")];I+=1
          if isnothing(X1) 
            X1,Y1=[x1,y1] 
          end

        # line to
        elseif cmd=="L" || cmd=="l"
          @assert(!isnothing(x1))
          x2,y2=[parse(Float64, jt) for jt in split_any(tokens[I]," ,")];I+=1
          push!(lines,[x1,y1, x2,y2])
          x1,y1=x2,y2

        # Close Path
        elseif cmd=="Z" || cmd=="z"
          @assert(!isnothing(X1))
          push!(lines,[x1,y1, X1,Y1])
        
        # Bezier Curves
        elseif cmd=="C" || cmd=="c"
          println("TODO cmd=[$(cmd)]")
          @assert(false) # unsupported, need to parse control points

        else
          println("Unsupported cmd=[$(cmd)]")
          @assert(false) # unsupported
        end
      end

      continue
    end

    # not supported
    assert(false)
  end

  return lines

end

# example
lines=parse_svg("test/book/svg/house-line.svg")

for (I,line) in enumerate(lines)
  println(I," ",line)
end


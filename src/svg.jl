
using Plasm
using EzXML
using DataStructures

# ///////////////////////////////////////////////////////////////
function cubicbezier2D(curvePts::Array{Array{Float64,1},1})
	b1(u) =       (1 - u)^3
	b2(u) = 3*u*  (1 - u)^2
	b3(u) = 3*u^2*(1 - u)
	b4(u) = u^3
	cntrlverts = hcat(curvePts...)
	x = cntrlverts[1,:]
	y = cntrlverts[2,:]
	Bx = u -> x[1]*b1(u) + x[2]*b2(u) + x[3]*b3(u) + x[4]*b4(u)
	By = u -> y[1]*b1(u) + y[2]*b2(u) + y[3]*b3(u) + y[4]*b4(u)
	return Bx,By
end

# ///////////////////////////////////////////////
function split_any(s, chars)::Vector{String}
	sep=string(chars[1])
	for other in chars[2:end]
		s=string(replace(s, string(other) => sep))
	end
	return Base.split(s,sep)
end

# ///////////////////////////////////////////////
function parse_svg_path(d::String)::Vector{String}

	if isempty(d)
		return []
	end

	if isletter(d[1])
		return [string(d[1]); parse_svg_path(d[2:end])]
	else
		I=1
		while I<=length(d) && !isletter(d[I]) I+=1 end
		return [strip(string(d[1:I-1])); parse_svg_path(d[I:end])]
	end
end


# ///////////////////////////////////////////////
function traverse_svg(cur, T::Matrix3d, points,hulls; bezier_show_control_points=false, bezier_npoints=10)

	function add_line(x1,y1,x2,y2)
		x1,y1,z1=T*[x1,y1,1.0];@assert(z1==1.0)
		x2,y2,z2=T*[x2,y2,1.0];@assert(z2==1.0)
		push!(points,[x1,y1])
		push!(points,[x2,y2])
		N=length(points)
		push!(hulls,[N-1,N])
	end

	# ___________________________________________
	if cur.name=="g"
		
		if haskey(cur,"transform")
			transform=cur["transform"]
			@assert(startswith(transform,"matrix("))
			@assert(endswith(transform,")"))
			transform=[parse(Float64, kt) for kt in split_any(transform[8:end-1],",")]
			T2=Matrix3d(
				transform[1],transform[3],transform[5],
				transform[2],transform[4],transform[6],
				0.0,0.0,1.0)
		else
			T2=Matrix3d()
		end

		for sub in eachelement(cur)
			traverse_svg(sub, T * T2, points, hulls, bezier_show_control_points=bezier_show_control_points, bezier_npoints=bezier_npoints)
		end
		return

	end

	# ___________________________________________
	if cur.name=="line"
		x1, y1, x2, y2=[parse(Float64, cur[k]) for k in ("x1","y1","x2","y2")]
		add_line(x1,y1, x2,y2)
		return
	end

	# ___________________________________________
	if cur.name=="rect"
		x1, y1, width, height=[parse(Float64, cur[k]) for k in ("x","y","width","height")]
		x2, y2=x1+width,y1+height
		add_line(x1,y1, x2,y1)
		add_line(x2,y1, x2,y2)
		add_line(x2,y2, x1,y2)
		add_line(x1,y2, x1,y1)
		return
	end

	# ___________________________________________
	if cur.name=="path"

		first_pos, current_pos=nothing,nothing
		tokens=parse_svg_path(cur["d"])

		Cursor=1
		function get_token()
			ret=tokens[Cursor]
			Cursor+=1
			return ret
		end
		
		while Cursor<=length(tokens)

			cmd=get_token()

			# Move To
			if cmd=="M" || cmd=="m"
				x,y=[parse(Float64, jt) for jt in split_any(get_token()," ,")]
				current_pos=[x,y]
				first_pos=current_pos

			# Line To
			elseif cmd=="L" || cmd=="l"
				@assert(!isnothing(current_pos))
				x,y=[parse(Float64, jt) for jt in split_any(get_token()," ,")]
				add_line(current_pos..., x,y)
				current_pos=[x,y]

			# Horizontal Line
			elseif cmd=="H" || cmd=="h"
				@assert(!isnothing(current_pos))
				x=parse(Float64, get_token())
				y=current_pos[2]
				add_line(current_pos..., x,y)
				current_pos=[x,y]

			# Vertical Line
		elseif cmd=="V" || cmd=="v"
			@assert(!isnothing(current_pos))
			x=current_pos[1]
			y=parse(Float64, get_token())
			add_line(current_pos..., x,y)
			current_pos=[x,y]

			# Close Path
			elseif cmd=="Z" || cmd=="z"
				@assert(!isnothing(first_pos) && !isnothing(current_pos))
				add_line(current_pos..., first_pos...)
				current_pos=first_pos
			
			# Bezier Curves
			elseif cmd=="C" || cmd=="c"
				@assert(!isnothing(current_pos))

				control_points=[current_pos]

				path_value=get_token()
				for jt in split_any(path_value," ")
					x,y=[parse(Float64, kt) for kt in split_any(jt,",")]
					push!(control_points,[x,y])
				end

				# todo, only cubic bezier supported right now
				if length(control_points)!=4
					println("SVG Bezier curve with ",length(control_points)," points. Ignoring since not supported")
					return
				end 

				Bx,By = cubicbezier2D(control_points)
				bezier_points = [ [Bx(u),By(u)] for u=0.0: (1.0/bezier_npoints) :1.0]

				for K=1:length(bezier_points)-1
					add_line(bezier_points[K]..., bezier_points[K+1]...)
				end

				if bezier_show_control_points
					add_line(control_points[1]..., control_points[2]...)
					add_line(control_points[2]..., control_points[3]...)
					add_line(control_points[3]..., control_points[4]...)
				end
				
				current_pos = control_points[end]

			else
				println("Unsupported cmd=[$(cmd)]")
				@assert(false) # unsupported
			end
		end

		return
	end


	# not supported, so ignoring
	println("SVG $(cur.name) not suppported so ignoring")

end


# ///////////////////////////////////////////////
function read_svg(filename::String; bezier_show_control_points=false, bezier_npoints=10)::Hpc

	points,hulls=Plasm.ConcretePointsNd(),Cells()
	doc = parsexml(read(filename,String))
	for it in eachelement(doc.root)
		traverse_svg(it, Matrix3d(), points, hulls, bezier_show_control_points=bezier_show_control_points, bezier_npoints=bezier_npoints)
	end


	return MKPOL(points,hulls)

end
export read_svg


"""
VIEW(read_svg("test/book/svg/house-line.svg"))
VIEW(read_svg("test/book/svg/Lar.svg"))
VIEW(read_svg("test/book/svg/curved.svg"))
VIEW(read_svg("test/book/svg/holes.svg"))
VIEW(read_svg("test/book/svg/new.svg"))
VIEW(read_svg("test/book/svg/tile.svg"))
"""

include("fenvs.jl")

# ///////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////
function Temple()

	function Column(r::Float64,h::Float64)
		basis = CUBOID([ 2*r*1.2, 2*r*1.2, h/12.0 ]) 
		trunk = CYLINDER([ r, (10.0/12.0)*h ])(12)
		capital = basis
		beam = S(1)(3)(capital) 
		return TOP([TOP([TOP([basis,trunk]),capital]),beam])
	end

	function Gable(radius::Float64,h::Float64,n::Int64)
		lastX = n*3*(2*radius*1.2)
		triangle = MKPOL(
			[[0,0],[lastX,0],[lastX/2.0,h/2]],
			[[1,2,3]],
			[[1]])
		return R([2,3])(PI/2)(POWER([triangle,QUOTE([radius*1.2])]))
	end

	col = Column(1.0, 12.0)

	function ColRow(N::Int64)
		v=[col for I in 1:N]
		return INSR(RIGHT)(v)
	end

	ColRowAndGable =  TOP([ColRow(4),Gable(1.0,12.0,4)])

	Temple = STRUCT(CAT([
		[ColRowAndGable, T(2)(6.0)], 
		DOUBLE_DIESIS(4)([ColRow(4),T(2)(6.0)]), 
		[ColRowAndGable] 
	]))
	
	Ground = EMBED(1)(BOX([1,2])(Temple))

	Xsizes = QUOTE( DOUBLE_DIESIS(14)([0.6,-1.2]) )
	Ysizes = QUOTE( AL([ -0.7, DOUBLE_DIESIS(5)([-1,5]) ]))
	Zsizes = QUOTE([ -13, 0.6 ])

	SecondaryBeams = POWER([POWER([Xsizes,Ysizes]),Zsizes])
	model= STRUCT([Temple, SecondaryBeams, Ground])
	VIEW(model)
end

# ////////////////////////////////////////////////
if abspath(PROGRAM_FILE) == @__FILE__
	Temple()
end
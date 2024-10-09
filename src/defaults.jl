const WHITE       = Point4d([1.0, 1.0, 1.0, 1.0])
const RED         = Point4d([1.0, 0.0, 0.0, 1.0])
const GREEN       = Point4d([0.0, 1.0, 0.0, 1.0])
const BLUE        = Point4d([0.0, 0.0, 1.0, 1.0])
const CYAN        = Point4d([0.0, 1.0, 1.0, 1.0])
const MAGENTA     = Point4d([1.0, 0.0, 1.0, 1.0])
const YELLOW      = Point4d([1.0, 1.0, 0.0, 1.0])
const ORANGE      = Point4d([1.0, 0.65, 1.0, 1.0])
const PURPLE      = Point4d([0.5, 0.0, 0.5, 1.0])
const BROWN       = Point4d([0.65, 0.16, 0.16, 1.0])
const GRAY        = Point4d([0.5, 0.5, 0.5, 1.0])
const BLACK       = Point4d([0.0, 0.0, 0.0, 1.0])
const DARK_GRAY   = Point4d(0.3,0.3,0.3,1.0)
const LIGHT_GRAY  = Point4d(0.8,0.8,0.8,1.0)
const TRANSPARENT = Point4d([0.0,0.0,0.0,0.0])

export WHITE,RED,GREEN,BLUE,CYAN,MAGENTA,YELLOW,ORANGE,PURPLE,BROWN,GRAY,BLACK, TRANSPARENT

function GetColorByName(name)
	return Dict{String,Point4d}(
		"white"=>WHITE, 
		"red"=>RED, 
		"green"=>GREEN, 
		"blue"=>BLUE,
		"cyan"=>CYAN, 
		"magenta"=>MAGENTA, 
		"yellow"=>YELLOW, 
		"orange"=>ORANGE,
		"purple"=>PURPLE, 
		"brown"=>BROWN, 
		"gray"=>GRAY, 
		"black"=>BLACK
	)
end
export GetColorByName


function RandomColor()::Point4d
	ret=Plasm.COLORS[rand(1:12)] # +0.5* Plasm.COLORS[rand(1:12)]
	#ret=1.0*[random_float(0.3,1.0) for I in 1:4]+0.0* Plasm.COLORS[rand(1:12)]
	ret[4] = 1.0
	return ret
end

DEFAULT_BACKGROUND_COLOR= Point4d([0.3,0.4,0.5,1.0]) 
DEFAULT_LAR_BACKGROUND_COLOR=WHITE

DEFAULT_USE_ORTHO=true
DEFAULT_SHOW_LINES=true
DEFAULT_SHOW_AXIS=true
DEFAULT_FOV=60.0
DEFAULT_LIGHTING_ENABLED=false

DEFAULT_POINT_SIZE   = 1
DEFAULT_LINE_WIDTH   = 3
DEFAULT_POINT_COLOR  = WHITE
DEFAULT_LINE_COLOR   = DARK_GRAY
DEFAULT_FACE_COLOR   = LIGHT_GRAY
DEFAULT_TEXT_SCALING=(0.250,0.375)

export DEFAULT_SHOW_AXIS,DEFAULT_BACKGROUND_COLOR,DEFAULT_LAR_BACKGROUND_COLOR, DEFAULT_SHOW_LINES,DEFAULT_USE_ORTHO, DEFAULT_FOV
export DEFAULT_POINT_SIZE,DEFAULT_POINT_COLOR
export DEFAULT_SHOW_LINES,DEFAULT_LINE_WIDTH,DEFAULT_LINE_COLOR
export DEFAULT_FACE_COLOR
export DEFAULT_LIGHTING_ENABLED
export DEFAULT_V_FONTSIZE,DEFAULT_EV_FONTSIZE,DEFAULT_FV_FONTSIZE
export DEFAULT_TEXT_SCALING

const MAYA01 = Point4d( 77/255.0, 202/255.0, 137/255.0, 1.0)
const MAYA02 = Point4d(130/255.0,  70/255.0,  88/255.0, 1.0)
const MAYA03 = Point4d(198/255.0, 180/255.0,  71/255.0, 1.0)
const MAYA04 = Point4d(133/255.0, 185/255.0,  98/255.0, 1.0)
const MAYA05 = Point4d(175/255.0, 115/255.0,  69/255.0, 1.0)
const MAYA06 = Point4d(227/255.0, 184/255.0, 128/255.0, 1.0)
const MAYA07 = Point4d(144/255.0, 121/255.0,  86/255.0, 1.0)
const MAYA08 = Point4d(203/255.0, 107/255.0, 113/255.0, 1.0)
const MAYA09 = Point4d(231/255.0, 157/255.0, 134/255.0, 1.0)
const MAYA10 = Point4d(146/255.0, 230/255.0, 136/255.0, 1.0)
const MAYA11 = Point4d( 79/255.0, 153/255.0, 125/255.0, 1.0)
const MAYA12 = Point4d(187/255.0, 154/255.0,  38/255.0, 1.0)

export MAYA01,MAYA02,MAYA03,MAYA04,MAYA05,MAYA06,MAYA07,MAYA08,MAYA09,MAYA10,MAYA11,MAYA12

const COLORS = [MAYA01,MAYA02,MAYA03,MAYA04,MAYA05,MAYA06,MAYA07,MAYA08,MAYA09,MAYA10,MAYA11,MAYA12]
export COLORS
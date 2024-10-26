using ModernGL
using GLFW


# /////////////////////////////////////////////////////////////////////
mutable struct GLVertexBuffer

	id::Int32
	vector::Vector{Float32}

  # constructor
  function GLVertexBuffer()
    return new(-1, Vector{Float32}())
  end

  # constructor
  function GLVertexBuffer(vector::Vector{Float32})
    return new(-1, copy(vector))
  end

end

export GLVertexBuffer

# /////////////////////////////////////////////////////////////////////
mutable struct GLVertexArray

	id::Int32

  # constructor
  function GLVertexArray()
    return new(-1)
  end

end
export GLVertexArray

# /////////////////////////////////////////////////////////////////////
mutable struct GLBatch
	primitive::UInt32
	T::Matrix4d
	vertex_array::GLVertexArray
	vertices::GLVertexBuffer
	normals::GLVertexBuffer
	colors::GLVertexBuffer
	point_size::Int
	line_width::Int
	point_color::Point4d
	line_color::Point4d
	face_color::Point4d
	enable_polygon_offset::Bool
end
export GLBatch


# ///////////////////////////////////////////////////////////////////////
function GLBatch(prim::UInt32=GL_POINTS)
	ret = GLBatch(
		prim,
		Matrix4d(), 
		GLVertexArray(),
		GLVertexBuffer(),
		GLVertexBuffer(),
		GLVertexBuffer(),
		1,
		1,
		Point4d(0.0, 0.0, 0.0, 1.0),
		Point4d(1.0, 1.0, 1.0, 1.0),
		Point4d(0.5, 0.5, 0.5, 1.0),
		false
	)
	return ret
end

# /////////////////////////////////////////////////////////////////////////////
mutable struct GLFWViewer  <: Viewer
	win::Any
	W::Int32
	H::Int32
	scalex::Float64
	scaley::Float64
	fov::Float64
	pos::Point3d
	dir::Point3d
	vup::Point3d
	zNear::Float64
	zFar::Float64
	walk_speed::Float64
	mouse_beginx::Float64
	mouse_beginy::Float64
	down_button::Int32
	batches::Vector{GLBatch}
	shaders::Dict{Any,Any}
	use_ortho::Bool
	exitNow::Bool
	show_lines::Bool
	background_color::PointNd
	title::String
	lighting_enabled::Bool
end
export GLFWViewer

# ///////////////////////////////////////////////////////////////////////
function GLFWViewer()
	ret=GLFWViewer(
		0,
		1024, 768,
		1.0, 1.0,
		DEFAULT_FOV,
		Point3d(0, 0, 1),
		Point3d(0, 0, -1),
		Point3d(1, 0, 0),
		0.01, 100.0, 0.1,
		0, 0, 0,
		Vector{GLBatch}(),
		Dict{Any,Any}(),
		DEFAULT_USE_ORTHO,
		false,  # exitNow
		DEFAULT_SHOW_LINES,
		DEFAULT_BACKGROUND_COLOR,
		"Plasm.jl",
		DEFAULT_LIGHTING_ENABLED
	)
	return ret
end



# /////////////////////////////////////////////////////////////////////
function glGenBuffer()
	id = GLuint[0]
	glGenBuffers(1, id)
	glCheckError("generating a buffer, array, or texture")
	id[]
end

# /////////////////////////////////////////////////////////////////////
function glGenVertexArray()
	id = GLuint[0]
	glGenVertexArrays(1, id)
	glCheckError("generating a buffer, array, or texture")
	id[]
end

# /////////////////////////////////////////////////////////////////////
function glCheckError(actionName="")
	message = glErrorMessage()
	if length(message) > 0
		if length(actionName) > 0
			error("Error ", actionName, ": ", message)
		else
			error("Error: ", message)
		end
	end
end

# /////////////////////////////////////////////////////////////////////
function glErrorMessage()
	err = glGetError()
	if err == GL_NO_ERROR
		return ""
	end
	if err == GL_INVALID_ENUM
		return "GL_INVALID_ENUM"
	end
	if err == GL_INVALID_VALUE
		return "GL_INVALID_VALUE"
	end
	if err == GL_INVALID_OPERATION
		return "GL_INVALID_OPERATION"
	end
	if err == GL_INVALID_FRAMEBUFFER_OPERATION
		return "GL_INVALID_FRAMEBUFFER_OPERATION"
	end
	if err == GL_OUT_OF_MEMORY
		return "GL_OUT_OF_MEMORY"
	end
	return "Unknown OpenGL error with error code"
end




# /////////////////////////////////////////////////////////////////////
function enableAttribute(location::Int32, buffer::GLVertexBuffer, num_components::Int64)
	if length(buffer.vector) == 00 || location < 0
		return
	end
	if buffer.id < 0
		buffer.id = glGenBuffer()
	end
	glBindBuffer(GL_ARRAY_BUFFER, buffer.id)
	glBufferData(GL_ARRAY_BUFFER, sizeof(buffer.vector), buffer.vector, GL_STATIC_DRAW)
	glVertexAttribPointer(location, num_components, GL_FLOAT, false, 0, C_NULL)
	glEnableVertexAttribArray(location)
	glBindBuffer(GL_ARRAY_BUFFER, 0)
end

# /////////////////////////////////////////////////////////////////////
function disableAttribute(location::Int32, buffer::GLVertexBuffer)
	if length(buffer.vector) == 00 || location < 0
		return
	end
	glDisableVertexAttribArray(location)
end


# /////////////////////////////////////////////////////////////////////
function enableVertexArray(array::GLVertexArray)
	if Sys.isapple()
		return
	end
	if array.id < 0
		array.id = glGenVertexArray()
	end
	glBindVertexArray(array.id)
end

# /////////////////////////////////////////////////////////////////////
function disableVertexArray(array::GLVertexArray)
	if Sys.isapple()
		return
	end
	glBindVertexArray(0)
end



# /////////////////////////////////////////////////////////////////////
mutable struct GLShader

	vertex_source
	frag_source

	program_id::Int32
	vertex_shader_id::Int32
	frag_shader_id::Int32

	# constructor
	function GLShader(vertex, fragment)
		return new(vertex, fragment, -1, -1, -1)
	end

end


# /////////////////////////////////////////////////////////////////////
function createShader(type, source)
	shader_id = glCreateShader(type)::GLuint
	glCheckError()
	glShaderSource(shader_id, 1, convert(Ptr{UInt8}, pointer([convert(Ptr{GLchar}, pointer(source))])), C_NULL)
	glCompileShader(shader_id)
	status = GLint[0]
	glGetShaderiv(shader_id, GL_COMPILE_STATUS, status)
	if status[1] == GL_FALSE
		maxlength = 8192
		buffer = zeros(GLchar, maxlength)
		sizei = GLsizei[0]
		glGetShaderInfoLog(shader_id, maxlength, sizei, buffer)
		len = sizei[]
		error_msg = unsafe_string(pointer(buffer), len)
		error("shader compilation failed\n", error_msg, "\nsource\n", source)
	end
	return shader_id
end


# /////////////////////////////////////////////////////////////////////
function enableProgram(shader)

	if (shader.program_id < 0)

		shader.program_id = glCreateProgram()
		glCheckError()

		shader.vertex_shader_id = createShader(GL_VERTEX_SHADER, shader.vertex_source)
		glAttachShader(shader.program_id, shader.vertex_shader_id)
		glCheckError()

		shader.frag_shader_id = createShader(GL_FRAGMENT_SHADER, shader.frag_source)
		glAttachShader(shader.program_id, shader.frag_shader_id)
		glCheckError()

		glLinkProgram(shader.program_id)
		glCheckError()
		status = GLint[0]
		glGetProgramiv(shader.program_id, GL_LINK_STATUS, status)
		if status[1] == GL_FALSE
			maxlength = 8192
			buffer = zeros(GLchar, maxlength)
			sizei = GLsizei[0]
			glGetProgramInfoLog(shader.program_id, maxlength, sizei, buffer)
			len = sizei[]
			error_msg = unsafe_string(pointer(buffer), len)
			error("Error linking program\n", error_msg)
		end
		glCheckError()
	end

	glUseProgram(shader.program_id)
end

# /////////////////////////////////////////////////////////////////////
function disableProgram(shader)
	glUseProgram(0)
end

# /////////////////////////////////////////////////////////////////////
phong_vert_source = """
	
	#define LIGHTING_ENABLED        arg(LIGHTING_ENABLED)
	#define COLOR_ATTRIBUTE_ENABLED arg(COLOR_ATTRIBUTE_ENABLED)
	
	uniform mat4 u_modelview_matrix;
	uniform mat4 u_projection_matrix;
	uniform vec4 u_color;
	
	attribute  vec4 a_position;
	
	#if LIGHTING_ENABLED
	attribute  vec3 a_normal;
	#endif
	
	#if COLOR_ATTRIBUTE_ENABLED
	attribute vec4 a_color;
	#endif
	
	#if LIGHTING_ENABLED
	uniform mat3 u_normal_matrix;
	uniform vec3 u_light_position;
	varying vec3 v_normal;
	varying vec3 v_light_dir;
	varying vec3 v_eye_vec;
	#endif
	
	#if COLOR_ATTRIBUTE_ENABLED
	varying vec4 v_color;
	#endif
	
	void main() 
	{
		vec4 eye_pos= u_modelview_matrix * a_position;
		
	#if LIGHTING_ENABLED	
		v_normal = u_normal_matrix * a_normal;
		vec3 vVertex = vec3(u_modelview_matrix * a_position);
		v_light_dir  = normalize(u_light_position - vVertex);
		v_eye_vec    = normalize(-vVertex);
	#endif	
	
	#if COLOR_ATTRIBUTE_ENABLED
		v_color=a_color;
	#endif
		
		gl_Position = u_projection_matrix * eye_pos;
	}
 """

# /////////////////////////////////////////////////////////////////////
phong_frag_source = """
	
	#define LIGHTING_ENABLED        arg(LIGHTING_ENABLED)
	#define COLOR_ATTRIBUTE_ENABLED arg(COLOR_ATTRIBUTE_ENABLED)

	uniform vec4 u_color;

	#if LIGHTING_ENABLED
	varying vec3 v_normal;
	varying vec3 v_light_dir;
	varying vec3 v_eye_vec;
	#endif

	#if COLOR_ATTRIBUTE_ENABLED
	varying vec4 v_color;
	#endif

	void main() 
	{
		vec4 frag_color=u_color; 
		
		#if LIGHTING_ENABLED
		vec3 N = normalize(v_normal   );
		vec3 L = normalize(v_light_dir);
		vec3 E = normalize(v_eye_vec  );

		vec4  u_material_ambient  = vec4(0.2,0.2,0.2,1.0);
		vec4  u_material_diffuse  = vec4(0.8,0.8,0.8,1.0) * u_color;
		vec4  u_material_specular = vec4(0.2,0.2,0.2,1.0) * u_color;
		float u_material_shininess=100.0;	
		
		if(gl_FrontFacing)
		{
			frag_color = u_material_ambient;
			float NdotL = abs(dot(N,L));
			if (NdotL>0.0)
				{
				vec3 R = reflect(-L, N);
				float NdotHV = abs(dot(R, E));
				frag_color += u_material_diffuse * NdotL;
				frag_color += u_material_specular * pow(NdotHV,u_material_shininess);
			}
		}
		else
		{
			frag_color = u_material_ambient;
			float NdotL = abs(dot(-N,L));
			if (NdotL>0.0);
			{
				vec3 R = reflect(-L, -N);
				float NdotHV=abs(dot(R, E));
				frag_color += u_material_diffuse * NdotL;
				frag_color += u_material_specular * pow(NdotHV,u_material_shininess);
			}
		}

	#if COLOR_ATTRIBUTE_ENABLED
	//frag_color = 0.5*v_color + frag_color*0.8;
	frag_color = v_color;
	#endif

	#elif COLOR_ATTRIBUTE_ENABLED
	frag_color = v_color;

	#endif

		gl_FragColor = frag_color;
	}
"""

# /////////////////////////////////////////////////////////////////////
function GLPhongShader(lighting_enabled, color_attribute_enabled)

	v = phong_vert_source
	f = phong_frag_source
	v = replace(v, "arg(LIGHTING_ENABLED)" => lighting_enabled ? "1" : "0")
	f = replace(f, "arg(LIGHTING_ENABLED)" => lighting_enabled ? "1" : "0")
	v = replace(v, "arg(COLOR_ATTRIBUTE_ENABLED)" => color_attribute_enabled ? "1" : "0")
	f = replace(f, "arg(COLOR_ATTRIBUTE_ENABLED)" => color_attribute_enabled ? "1" : "0")
	return GLShader(v, f)
end




# ///////////////////////////////////////////////////////////////////////
function getModelview(viewer::GLFWViewer)
	return lookAtMatrix(viewer.pos, viewer.pos + viewer.dir, viewer.vup)
end

function getProjection(viewer::GLFWViewer)
	ratio = viewer.W / float(viewer.H)
	if viewer.use_ortho
		# euristic that seem to work well
		Z = viewer.zNear + 0.5 * (viewer.zFar - viewer.zNear)
		right = Z * tan(deg2rad(viewer.fov / 2.0))
		left = -right
		top =    -0.5 * (right - left) / ratio
		bottom = +0.5 * (right - left) / ratio
		return orthoMatrix(left, right, top, bottom, viewer.zNear, viewer.zFar)
	else
		return perspectiveMatrix(viewer.fov, ratio, viewer.zNear, viewer.zFar)
	end

end

# ///////////////////////////////////////////////////////////////////////
function getShader(viewer::GLFWViewer, lighting_enabled, color_attribute_enabled)

	lighting_enabled = lighting_enabled && viewer.lighting_enabled

	key = (lighting_enabled, color_attribute_enabled)

	if haskey(viewer.shaders, key)
		return viewer.shaders[key]
	end

	ret = GLPhongShader(lighting_enabled, color_attribute_enabled)
	viewer.shaders[key] = ret
	return ret
end

# ///////////////////////////////////////////////////////////////////////
function glRender(viewer::GLFWViewer)

	glDeleteNow()

	glEnable(GL_BLEND)
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

	glEnable(GL_DEPTH_TEST)
	glDepthFunc(GL_LEQUAL)
	glDisable(GL_CULL_FACE)
	glClearDepth(1.0)
	glClearColor(viewer.background_color[1], viewer.background_color[2], viewer.background_color[3], 0.00)
	glPolygonOffset(-1.0, -1.0)

	glViewport(0, 0, viewer.W, viewer.H)
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

	PROJECTION = getProjection(viewer)
	MODELVIEW = getModelview(viewer)

	lightpos = MODELVIEW * Point4d(viewer.pos[1], viewer.pos[2], viewer.pos[3], 1.0)

	for batch in viewer.batches

		# show points
		if batch.primitive == GL_POINTS && (batch.point_color[4] > 0.0)
			glPointSize(batch.point_size)
			render_batch(viewer, batch, batch.point_color, PROJECTION, MODELVIEW, lightpos)
			glPointSize(1)

		# show lines
		elseif (batch.primitive == GL_LINES || batch.primitive == GL_LINE_STRIP || batch.primitive == GL_LINE_LOOP) && (batch.line_color[4] > 0.0)
			glLineWidth(batch.line_width)
			render_batch(viewer, batch, batch.line_color, PROJECTION, MODELVIEW, lightpos)
			glLineWidth(1)

		# show triangles
		elseif (batch.primitive == GL_TRIANGLES) && (batch.face_color[4] > 0.0)
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
			if !isnothing(batch.enable_polygon_offset)
				glEnable(GL_POLYGON_OFFSET_FILL)
				glPolygonOffset(1.0, 1.0)
				render_batch(viewer, batch, batch.face_color, PROJECTION, MODELVIEW, lightpos)
				glDisable(GL_POLYGON_OFFSET_FILL)
			else
				render_batch(viewer, batch, batch.face_color, PROJECTION, MODELVIEW, lightpos)
			end
		end
	end

	glCheckError()

end

# //////////////////////////////////////////////////////////////////////////////////////////////
function render_batch(viewer::GLFWViewer, batch::GLBatch, color, PROJECTION, MODELVIEW, lightpos)

	lighting_enabled = length(batch.normals.vector) > 0 && viewer.lighting_enabled
	color_attribute_enabled = length(batch.colors.vector) > 0

	shader = getShader(viewer, lighting_enabled, color_attribute_enabled)

	enableProgram(shader)

	projection = PROJECTION
	modelview = MODELVIEW * batch.T
	normal_matrix = dropW(transpose(inv(modelview)))

	glUniformMatrix4fv(glGetUniformLocation(shader.program_id, "u_modelview_matrix"), 1, GL_TRUE, flatten(modelview))
	glUniformMatrix4fv(glGetUniformLocation(shader.program_id, "u_projection_matrix"), 1, GL_TRUE, flatten(projection))
	glUniformMatrix3fv(glGetUniformLocation(shader.program_id, "u_normal_matrix"), 1, GL_TRUE, flatten(normal_matrix))

	u_light_position = glGetUniformLocation(shader.program_id, "u_light_position")
	if u_light_position >= 0
		glUniform3f(u_light_position, lightpos[1] / lightpos[4], lightpos[2] / lightpos[4], lightpos[3] / lightpos[4])
	end

	u_color = glGetUniformLocation(shader.program_id, "u_color")
	# println(u_color," ", color, " lighting_enabled ",lighting_enabled, " color_attribute_enabled ",color_attribute_enabled)
	if u_color >= 0 && !isnothing(color)
		glUniform4f(u_color, color[1], color[2], color[3], color[4])
	end

	enableVertexArray(batch.vertex_array)

	a_position = glGetAttribLocation(shader.program_id, "a_position")
	a_normal = glGetAttribLocation(shader.program_id, "a_normal")
	a_color = glGetAttribLocation(shader.program_id, "a_color")

	enableAttribute(a_position, batch.vertices, 3)
	enableAttribute(a_normal, batch.normals, 3)
	enableAttribute(a_color, batch.colors, 4)

	@assert length(batch.vertices.vector) % 3 == 0
	glDrawArrays(batch.primitive, 0, Int64(length(batch.vertices.vector) / 3))

	disableAttribute(a_position, batch.vertices)
	disableAttribute(a_normal, batch.normals)
	disableAttribute(a_color, batch.colors)
	disableVertexArray(batch.vertex_array)
	disableProgram(shader)

	glDepthMask(true)
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
	glDisable(GL_POLYGON_OFFSET_LINE)

end

# ///////////////////////////////////////////////////////////////////////
function redisplay(viewer::GLFWViewer)
	# nothing to do
end

# ///////////////////////////////////////////////////////////////////////
function handleResizeEvent(viewer)
	size = GLFW.GetWindowSize(viewer.win)
	viewer.W = size[1] * viewer.scalex
	viewer.H = size[2] * viewer.scaley
	redisplay(viewer)
end

# ///////////////////////////////////////////////////////////////////////
function handleMouseButtonEvent(viewer::GLFWViewer, button, action, mods)

	button = Dict(GLFW.MOUSE_BUTTON_1 => 1, GLFW.MOUSE_BUTTON_2 => 3, GLFW.MOUSE_BUTTON_3 => 2)[button]

	if action == GLFW.PRESS && viewer.down_button == 0
		viewer.down_button = button
		redisplay(viewer)
		return
	end

	if action == GLFW.RELEASE && button == viewer.down_button
		viewer.down_button = 0
		redisplay(viewer)
		return
	end
end

# ///////////////////////////////////////////////////////////////////////
function handleMouseMoveEvent(viewer::GLFWViewer, x::Float64, y::Float64)

	x = x * viewer.scalex
	y = y * viewer.scaley

	button = viewer.down_button

	if (button == 0)
		viewer.mouse_beginx = x
		viewer.mouse_beginy = y
		return
	end

	deltax = float(x - viewer.mouse_beginx)
	deltay = float(viewer.mouse_beginy - y)
	W = viewer.W
	H = viewer.H

	modelview = getModelview(viewer)

	if button == 1
		screen_center = Point3d(W / 2.0, H / 2.0, 0.0)
		a = (Point3d((float)(viewer.mouse_beginx - screen_center[1]), (float)(H - viewer.mouse_beginy - screen_center[2]), 0)) * (1.0 / min(W, H))
		b = (Point3d((float)(x - screen_center[1]), (float)(H - y - screen_center[2]), 0)) * (1.0 / min(W, H))
		a[3] = 2.0^(-0.5 * norm(a))
		b[3] = 2.0^(-0.5 * norm(b))
		a = normalized(a)
		b = normalized(b)
		axis = normalized(cross(a, b))
		angle = acos(dot(a, b))

		#vt=Point3d(modelview[1,4],modelview[2,4],modelview[3,4])
		#modelview=translateMatrix(vt) * convertToMatrix(convertToQuaternion(modelview))

		q = Quaternion(axis, angle) * convertToQuaternion(modelview)
		vt = Point3d(modelview[1, 4], modelview[2, 4], modelview[3, 4])
		modelview = translateMatrix(vt) * convertToMatrix(q)

	elseif button == 3
		vt = Point3d(deltax * viewer.walk_speed, deltay * viewer.walk_speed, 0.0)
		modelview = translateMatrix(vt) * modelview
	end

	viewer.pos, viewer.dir, viewer.vup = getLookAt(modelview)

	viewer.mouse_beginx = x
	viewer.mouse_beginy = y
	redisplay(viewer)
end

# ///////////////////////////////////////////////////////////////////////
function handleMouseWheelEvent(viewer::GLFWViewer, delta)
	if viewer.use_ortho
		viewer.fov*=(delta<0) ? 1.1 : 1.0/1.1;
	else
		viewer.pos = viewer.pos + viewer.dir * ((delta >= 0 ? 10.0 : -10.0) * viewer.walk_speed)
	end
	redisplay(viewer)
end

# ///////////////////////////////////////////////////////////////////////
function handleKeyPressEvent(viewer::GLFWViewer, key, scancode, action, mods)

	if action != GLFW.PRESS && action != GLFW.REPEAT
		return
	end

	if key == GLFW.KEY_ESCAPE
		viewer.exitNow = true
		return
	end

	if (key == GLFW.KEY_KP_ADD)
		viewer.walk_speed *= 0.95
		return
	end

	if (key == GLFW.KEY_KP_SUBTRACT)
		viewer.walk_speed *= (1.0 / 0.95)
		return
	end

	if (key == GLFW.KEY_W)
		viewer.pos = viewer.pos + viewer.dir * viewer.walk_speed
		redisplay(viewer)
		return
	end

	if (key == GLFW.KEY_S)
		viewer.pos = viewer.pos - viewer.dir * viewer.walk_speed
		redisplay(viewer)
		return
	end

	if (key == GLFW.KEY_O)
		viewer.use_ortho = !viewer.use_ortho
		#println("use_ortho ",viewer.use_ortho)
		redisplay(viewer)
		return
	end


	if (key == GLFW.KEY_UP)
		viewer.pos = viewer.pos + viewer.vup * viewer.walk_speed
		redisplay(viewer)
		return
	end

	if (key == GLFW.KEY_DOWN)
		viewer.pos = viewer.pos - viewer.vup * viewer.walk_speed
		redisplay(viewer)
		return
	end

	if (key == GLFW.KEY_LEFT || key == GLFW.KEY_A)
		right = normalized(cross(viewer.dir, viewer.vup))
		viewer.pos = viewer.pos - right * viewer.walk_speed
		redisplay(viewer)
		return
	end

	if (key == GLFW.KEY_RIGHT || key == GLFW.KEY_D)
		right = normalized(cross(viewer.dir, viewer.vup))
		viewer.pos = viewer.pos + right * viewer.walk_speed
		redisplay(viewer)
		return
	end

	if (key == GLFW.KEY_L)
		viewer.show_lines = !viewer.show_lines
		redisplay(viewer)
		return
	end

end

# ///////////////////////////////////////////////////////////////////////
function compute_bbox(batch::GLBatch)::Box3d
	ret = invalidBox()
	for (x,y,z) in split(batch.vertices.vector,3)
		addPoint(ret, Point3d(x,y,z))
	end
	return ret
end

# ///////////////////////////////////////////////////////////////////////
function compute_bbox(batches::Vector{GLBatch})::Box3d
	ret = invalidBox()
	for batch in batches
		for (x,y,z) in split(batch.vertices.vector,3)
			addPoint(ret, Point3d(x,y,z))
		end
	end
	return ret
end

# ///////////////////////////////////////////////////////////////////////
function render_points(viewer::GLFWViewer, points::Vector{Float32}; colors=nothing, point_size=DEFAULT_POINT_SIZE, point_color=DEFAULT_POINT_COLOR)
	if length(points)==0 return end
	batch = GLBatch(GL_POINTS)
	push!(viewer.batches,batch)
	batch.point_size=point_size
	point_color=point_color
	batch.vertices = GLVertexBuffer(points)
	if !isnothing(colors)
		batch.colors = GLVertexBuffer(colors)
	end
end

# ///////////////////////////////////////////////////////////////////////
function render_lines(viewer::GLFWViewer, points::Vector{Float32}; colors=nothing, line_width=DEFAULT_LINE_WIDTH, line_color=DEFAULT_LINE_COLOR)
	if length(points)==0 return end
	batch = GLBatch(GL_LINES)
	push!(viewer.batches, batch)
	batch.line_width=line_width 
	batch.line_color=line_color
	batch.vertices = GLVertexBuffer(points)
	if !isnothing(colors)
		batch.colors = GLVertexBuffer(colors)
	end
end

# ///////////////////////////////////////////////////////////////////////
function render_triangles(viewer::GLFWViewer, points::Vector{Float32}; normals=nothing, colors=nothing, face_color=DEFAULT_FACE_COLOR, enable_polygon_offset=false)
	if length(points)==0 return end
	batch = GLBatch(GL_TRIANGLES)
	push!(viewer.batches,batch)
	batch.enable_polygon_offset=enable_polygon_offset
	batch.face_color=face_color
	batch.vertices=GLVertexBuffer(points)
	batch.normals=GLVertexBuffer(!isnothing(normals) ? normals : compute_triangles_normals(points))
	if !isnothing(colors)
		batch.colors = GLVertexBuffer(colors)
		
	end
end


# /////////////////////////////////////////////////////////////////////
__release_gpu_resources__ = []

function glDeleteLater(fun::Function)
	global __release_gpu_resources__
	append!(__release_gpu_resources__, [fun])
end

function glDeleteNow()
	global __release_gpu_resources__
	v,__release_gpu_resources__=__release_gpu_resources__,[]
	for fun in __release_gpu_resources__
		fun()
	end
end 

function releaseGpuResources(array::GLVertexArray)
       global __release_gpu_resources__
       if array.id >= 0
               id = array.id
               array.id = -1
               glDeleteLater(function () glDeleteVertexArrays(1, [id])end)
       end
end

function releaseGpuResources(buffer::GLVertexBuffer)
	global __release_gpu_resources__
	if buffer.id >= 0
		id = buffer.id
		buffer.id = -1
		glDeleteLater(function () glDeleteBuffers(1, [id])end)
	end
end
function releaseGpuResources(batch::GLBatch)
	releaseGpuResources(batch.vertex_array)
	releaseGpuResources(batch.vertices)
	releaseGpuResources(batch.normals)
	releaseGpuResources(batch.colors)
end

function releaseGpuResources(shader::GLShader)
	global __release_gpu_resources__

	if shader.vertex_shader_id >= 0
					id = shader.vertex_shader_id
					shader.vertex_shader_id = -1
					glDeleteLater(function () glDeleteShader(id)end)
	end

	if shader.frag_shader_id >= 0
					id = shader.frag_shader_id
					shader.frag_shader_id = -1
					glDeleteLater(function ()glDeleteShader(id)end)

	end

	if shader.program_id >= 0
					id = shader.program_id
					shader.program_id = -1
					glDeleteLater(function ()glDeleteProgram(id)end)
	end

end


# ///////////////////////////////////////////////////////////////////////
function run_viewer_in_current_thread(viewer::GLFWViewer; properties::Properties=Properties())

	# there should NOT be any repetition
	@assert(length(Set(viewer.batches))==length(viewer.batches))

	# compute bbox before adding the axis
	begin
		bbox=compute_bbox(viewer.batches)
		box_size = bbox.p2 - bbox.p1
		max_size = maximum([box_size[1], box_size[2], box_size[3]])
		box_center = center(bbox)
		use_ortho=(box_size[3]==0.0)
	end

	show_axis = get(properties, "show_axis", true)
	if show_axis
		render_axis(viewer, Point3d(0, 0, 0), Point3d(1.05, 1.05, 1.05))
	end

	# reduce to the case -1,+1 (this is needed for meshcat)
	if true
		T =  scaleMatrix(Point3d(2.0/max_size, 2.0/max_size, 2.0/max_size)) * translateMatrix(-box_center)
		for (B,batch) in enumerate(viewer.batches)

			batch.T= T * batch.T # prepend transformation

			# batch.vertices.vector=transform_points(T, batch.vertices.vector)
			# note: scaling/translation will not change the normals so I will be fine
			# batch.normals.vector=transform_normals(T,batch.normals.vector)

		end
		
	end

	viewer.background_color =            get(properties, "background_color", viewer.background_color)
	viewer.title            =            get(properties, "title", viewer.title)
	viewer.use_ortho        =            get(properties, "use_ortho", use_ortho)
	viewer.show_lines       =            get(properties, "show_lines", viewer.show_lines)
	viewer.fov              =            get(properties, "fov", viewer.fov)

	viewer.pos              =            use_ortho ? Point3d(0, 0, 1) : Point3d(5, 5, 5) 
	viewer.dir              = normalized(normalized(Point3d(0,0,0) - viewer.pos))
	viewer.vup              = normalized(use_ortho ? Point3d(0, 1, 0) : Point3d(0, 0, 1))
	viewer.zNear            =            0.01
	viewer.zFar             =            20.0
	viewer.walk_speed       =            0.01 

	redisplay(viewer)

	ret_code = GLFW.Init()
	# println("GLFW init returned ",ret_code)

	win = GLFW.CreateWindow(viewer.W, viewer.H, viewer.title)
	viewer.win = win
	viewer.exitNow = false
	GLFW.MakeContextCurrent(win)

	# problem of retina
	window_size = GLFW.GetWindowSize(viewer.win)
	framebuffer_size = GLFW.GetFramebufferSize(viewer.win)
	viewer.scalex = framebuffer_size[1] / Float64(window_size[1])
	viewer.scaley = framebuffer_size[2] / Float64(window_size[2])

	GLFW.SetWindowSizeCallback(win, function (win::GLFW.Window, width::Integer, height::Integer)
		handleResizeEvent(viewer)
	end)
	GLFW.SetKeyCallback(win, function (win::GLFW.Window, key, scancode, action, mods)
		handleKeyPressEvent(viewer, key, scancode, action, mods)
	end)
	GLFW.SetCursorPosCallback(win, function (win::GLFW.Window, x, y)
		handleMouseMoveEvent(viewer, x, y)
	end)
	GLFW.SetMouseButtonCallback(win, function (win::GLFW.Window, button, action, mods)
		handleMouseButtonEvent(viewer, button, action, mods)
	end)
	GLFW.SetScrollCallback(win, function (win::GLFW.Window, dx, dy)
		handleMouseWheelEvent(viewer, dy)
	end)

	handleResizeEvent(viewer::GLFWViewer)

	while !GLFW.WindowShouldClose(win)
		if viewer.exitNow SetWindowShouldClose(win,true) end
		glRender(viewer)
		GLFW.SwapBuffers(win)
		GLFW.PollEvents()
		sleep(0.001)
	end

	for batch in viewer.batches
		releaseGpuResources(batch)
	end
	
	# keeping shaders
	# for (key, shader) in viewer.shaders
	#	 releaseGpuResources(shader)
	#end

	# GLFW.HideWindow(win)
	GLFW.DestroyWindow(win)
	
	# println("Viewer destroyed")
end


# ///////////////////////////////////////////////////////////////////////
function run_viewer(viewer::GLFWViewer; properties::Properties=Properties(), use_thread=true)
	#println("Creating Window use_thread=",use_thread)
	#if use_thread
	#	Threads.@spawn run_viewer_in_current_thread(viewer, properties=properties)
	#else
		run_viewer_in_current_thread(viewer, properties=properties)
		GLFW.Terminate()
	#end
end
export run_viewer
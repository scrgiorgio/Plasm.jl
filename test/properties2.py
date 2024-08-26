
using Plasm

# ////////////////////////////////////////////////////////
# geometry (per-batch) properties 
# ///////////////////////////////////////////////////////
obj=PROPERTIES(
  CUBOID([1,1,1]),
  Properties(
    "point_size"        => DEFAULT_POINT_SIZE,
    "point_color"       => DEFAULT_POINT_COLOR,
    "line_width"        => DEFAULT_LINE_WIDTH,
    "line_color"        => DEFAULT_LINE_COLOR, 
    "face_color"        => DEFAULT_FACE_COLOR,

    # used by certain algorithms to set/retrieve node type
    "node_type"         =>"whatever-you-need",
  ))

# ////////////////////////////////////////////////////////
# viewer (global) properties
# ////////////////////////////////////////////////////////
viewer_properties=Properties(

  "title"            => "My viewer example",

  "background_color" => DEFAULT_BACKGROUND_COLOR,

  "show_axis"        => DEFAULT_SHOW_AXIS,
  "show_lines"       => DEFAULT_SHOW_LINES,
  "lighting_enabled" => DEFAULT_LIGHTING_ENABLED,

  "use_ortho"        => DEFAULT_USE_ORTHO,

  # example viewer at Z=2 looking downward at origin (X is vup) 
  "pos"              => Point3d(0,0,+2),
  "dir"              => Point3d(0,0,-1),
  "vup"              => Point3d(1,0,0),
  "fov"              => DEFAULT_FOV,

  # opengl camera parameters
  "znear"            => 0.1,
  "zfar"             => 100.0,
  "walk_speed"       => 10.0,

  # LAR VIEWCOMPLEX  specific
  "text_v_color"      => DEFAULT_TEXT_V_COLOR,
  "text_ev_color"     => DEFAULT_TEXT_EV_COLOR,
  "text_fv_color"     => DEFAULT_TEXT_FV_COLOR,

  "v_fontsize"        => DEFAULT_V_FONTSIZE,
  "ev_fontsize"       => DEFAULT_EV_FONTSIZE,
  "fv_fontsize"       => DEFAULT_FV_FONTSIZE,
)

VIEW(obj, viewer_properties)


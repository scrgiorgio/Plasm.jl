
using XML
using LinearAlgebra
using Printf

export parse_svg_file, convert_svg_to_paths


# ////////////////////////////////////////////////////////////////////////
function parse_color(color_str::String)::Vector{Float64}
    rgba = [0.0, 0.0, 0.0, 1.0]
    if isempty(color_str) || lowercase(color_str) == "none"
        return [0.0, 0.0, 0.0, 0.0]  # Transparent
    end
    if startswith(color_str, "#")
        return parse_hex_color(color_str)
    end
    if startswith(lowercase(color_str), "rgb")
        return parse_rgb_color(color_str)
    end
    named_color = get_named_color(color_str)
    if !isnothing(named_color)
        return named_color
    end
    return rgba
end


# ////////////////////////////////////////////////////////////////////////
function parse_hex_color(hex_str::String)::Vector{Float64}
    hex = replace(hex_str, "#" => "")
    if length(hex) == 3
        r = parse(Int, hex[1] * hex[1], base=16) / 255.0
        g = parse(Int, hex[2] * hex[2], base=16) / 255.0
        b = parse(Int, hex[3] * hex[3], base=16) / 255.0
        return [r, g, b, 1.0]
    end
    if length(hex) == 4
        r = parse(Int, hex[1] * hex[1], base=16) / 255.0
        g = parse(Int, hex[2] * hex[2], base=16) / 255.0
        b = parse(Int, hex[3] * hex[3], base=16) / 255.0
        a = parse(Int, hex[4] * hex[4], base=16) / 255.0
        return [r, g, b, a]
    end
    if length(hex) == 6
        r = parse(Int, hex[1:2], base=16) / 255.0
        g = parse(Int, hex[3:4], base=16) / 255.0
        b = parse(Int, hex[5:6], base=16) / 255.0
        return [r, g, b, 1.0]
    end
    if length(hex) == 8
        r = parse(Int, hex[1:2], base=16) / 255.0
        g = parse(Int, hex[3:4], base=16) / 255.0
        b = parse(Int, hex[5:6], base=16) / 255.0
        a = parse(Int, hex[7:8], base=16) / 255.0
        return [r, g, b, a]
    end
    return [0.0, 0.0, 0.0, 1.0]
end


# ////////////////////////////////////////////////////////////////////////
function parse_rgb_color(rgb_str::String)::Vector{Float64}
    is_rgba = startswith(lowercase(rgb_str), "rgba")
    values_str = match(r"\((.*)\)", rgb_str).captures[1]
    values = Base.split(values_str, r"[,\s]+")
    r = parse_color_component(String(strip(values[1])))
    g = parse_color_component(String(strip(values[2])))
    b = parse_color_component(String(strip(values[3])))
    a = 1.0
    if is_rgba && length(values) >= 4
        a = parse(Float64, strip(values[4]))
    end
    
    return [r, g, b, a]
end


# ////////////////////////////////////////////////////////////////////////
function parse_color_component(component_str::AbstractString)::Float64
    component_str = strip(component_str)
    if endswith(component_str, "%")
        percentage = parse(Float64, replace(component_str, "%" => ""))
        return clamp(percentage / 100.0, 0.0, 1.0)
    end
    value = parse(Float64, component_str)
    if value <= 1.0 && !contains(component_str, ".")
        return clamp(value / 255.0, 0.0, 1.0)
    else
        return clamp(value, 0.0, 1.0)
    end
end


# ////////////////////////////////////////////////////////////////////////
function get_named_color(name::String)::Union{Vector{Float64}, Nothing}
    named_colors = Dict{String, Vector{Float64}}(
        "black" => [0.0, 0.0, 0.0, 1.0],
        "silver" => [0.75, 0.75, 0.75, 1.0],
        "gray" => [0.5, 0.5, 0.5, 1.0],
        "white" => [1.0, 1.0, 1.0, 1.0],
        "maroon" => [0.5, 0.0, 0.0, 1.0],
        "red" => [1.0, 0.0, 0.0, 1.0],
        "purple" => [0.5, 0.0, 0.5, 1.0],
        "fuchsia" => [1.0, 0.0, 1.0, 1.0],
        "green" => [0.0, 0.5, 0.0, 1.0],
        "lime" => [0.0, 1.0, 0.0, 1.0],
        "olive" => [0.5, 0.5, 0.0, 1.0],
        "yellow" => [1.0, 1.0, 0.0, 1.0],
        "navy" => [0.0, 0.0, 0.5, 1.0],
        "blue" => [0.0, 0.0, 1.0, 1.0],
        "teal" => [0.0, 0.5, 0.5, 1.0],
        "aqua" => [0.0, 1.0, 1.0, 1.0],
        "orange" => [1.0, 0.65, 0.0, 1.0],
        "aliceblue" => [0.94, 0.97, 1.0, 1.0],
        "antiquewhite" => [0.98, 0.92, 0.84, 1.0],
        "aquamarine" => [0.5, 1.0, 0.83, 1.0],
        "azure" => [0.94, 1.0, 1.0, 1.0],
        "beige" => [0.96, 0.96, 0.86, 1.0],
        "bisque" => [1.0, 0.89, 0.77, 1.0],
        "blanchedalmond" => [1.0, 0.92, 0.8, 1.0],
        "blueviolet" => [0.54, 0.17, 0.89, 1.0],
        "brown" => [0.65, 0.16, 0.16, 1.0],
        "burlywood" => [0.87, 0.72, 0.53, 1.0],
        "cadetblue" => [0.37, 0.62, 0.63, 1.0],
        "chartreuse" => [0.5, 1.0, 0.0, 1.0],
        "chocolate" => [0.82, 0.41, 0.12, 1.0],
        "coral" => [1.0, 0.5, 0.31, 1.0],
        "cornflowerblue" => [0.39, 0.58, 0.93, 1.0],
        "cornsilk" => [1.0, 0.97, 0.86, 1.0],
        "crimson" => [0.86, 0.08, 0.24, 1.0],
        "cyan" => [0.0, 1.0, 1.0, 1.0],
        "darkblue" => [0.0, 0.0, 0.55, 1.0],
        "darkcyan" => [0.0, 0.55, 0.55, 1.0],
        "darkgoldenrod" => [0.72, 0.53, 0.04, 1.0],
        "darkgray" => [0.66, 0.66, 0.66, 1.0],
        "darkgreen" => [0.0, 0.39, 0.0, 1.0],
        "darkgrey" => [0.66, 0.66, 0.66, 1.0],
        "darkkhaki" => [0.74, 0.72, 0.42, 1.0],
        "darkmagenta" => [0.55, 0.0, 0.55, 1.0],
        "darkolivegreen" => [0.33, 0.42, 0.18, 1.0],
        "darkorange" => [1.0, 0.55, 0.0, 1.0],
        "darkorchid" => [0.6, 0.2, 0.8, 1.0],
        "darkred" => [0.55, 0.0, 0.0, 1.0],
        "darksalmon" => [0.91, 0.59, 0.48, 1.0],
        "darkseagreen" => [0.56, 0.74, 0.56, 1.0],
        "darkslateblue" => [0.28, 0.24, 0.55, 1.0],
        "darkslategray" => [0.18, 0.31, 0.31, 1.0],
        "darkslategrey" => [0.18, 0.31, 0.31, 1.0],
        "darkturquoise" => [0.0, 0.81, 0.82, 1.0],
        "darkviolet" => [0.58, 0.0, 0.83, 1.0],
        "deeppink" => [1.0, 0.08, 0.58, 1.0],
        "deepskyblue" => [0.0, 0.75, 1.0, 1.0],
        "dimgray" => [0.41, 0.41, 0.41, 1.0],
        "dimgrey" => [0.41, 0.41, 0.41, 1.0],
        "dodgerblue" => [0.12, 0.56, 1.0, 1.0],
        "firebrick" => [0.7, 0.13, 0.13, 1.0],
        "floralwhite" => [1.0, 0.98, 0.94, 1.0],
        "forestgreen" => [0.13, 0.55, 0.13, 1.0],
        "gainsboro" => [0.86, 0.86, 0.86, 1.0],
        "ghostwhite" => [0.97, 0.97, 1.0, 1.0],
        "gold" => [1.0, 0.84, 0.0, 1.0],
        "goldenrod" => [0.85, 0.65, 0.13, 1.0],
        "greenyellow" => [0.68, 1.0, 0.18, 1.0],
        "grey" => [0.5, 0.5, 0.5, 1.0],
        "honeydew" => [0.94, 1.0, 0.94, 1.0],
        "hotpink" => [1.0, 0.41, 0.71, 1.0],
        "indianred" => [0.8, 0.36, 0.36, 1.0],
        "indigo" => [0.29, 0.0, 0.51, 1.0],
        "ivory" => [1.0, 1.0, 0.94, 1.0],
        "khaki" => [0.94, 0.9, 0.55, 1.0],
        "lavender" => [0.9, 0.9, 0.98, 1.0],
        "lavenderblush" => [1.0, 0.94, 0.96, 1.0],
        "lawngreen" => [0.49, 0.99, 0.0, 1.0],
        "lemonchiffon" => [1.0, 0.98, 0.8, 1.0],
        "lightblue" => [0.68, 0.85, 0.9, 1.0],
        "lightcoral" => [0.94, 0.5, 0.5, 1.0],
        "lightcyan" => [0.88, 1.0, 1.0, 1.0],
        "lightgoldenrodyellow" => [0.98, 0.98, 0.82, 1.0],
        "lightgray" => [0.83, 0.83, 0.83, 1.0],
        "lightgreen" => [0.56, 0.93, 0.56, 1.0],
        "lightgrey" => [0.83, 0.83, 0.83, 1.0],
        "lightpink" => [1.0, 0.71, 0.76, 1.0],
        "lightsalmon" => [1.0, 0.63, 0.48, 1.0],
        "lightseagreen" => [0.13, 0.7, 0.67, 1.0],
        "lightskyblue" => [0.53, 0.81, 0.98, 1.0],
        "lightslategray" => [0.47, 0.53, 0.6, 1.0],
        "lightslategrey" => [0.47, 0.53, 0.6, 1.0],
        "lightsteelblue" => [0.69, 0.77, 0.87, 1.0],
        "lightyellow" => [1.0, 1.0, 0.88, 1.0],
        "limegreen" => [0.2, 0.8, 0.2, 1.0],
        "linen" => [0.98, 0.94, 0.9, 1.0],
        "magenta" => [1.0, 0.0, 1.0, 1.0],
        "mediumaquamarine" => [0.4, 0.8, 0.67, 1.0],
        "mediumblue" => [0.0, 0.0, 0.8, 1.0],
        "mediumorchid" => [0.73, 0.33, 0.83, 1.0],
        "mediumpurple" => [0.58, 0.44, 0.86, 1.0],
        "mediumseagreen" => [0.24, 0.7, 0.44, 1.0],
        "mediumslateblue" => [0.48, 0.41, 0.93, 1.0],
        "mediumspringgreen" => [0.0, 0.98, 0.6, 1.0],
        "mediumturquoise" => [0.28, 0.82, 0.8, 1.0],
        "mediumvioletred" => [0.78, 0.08, 0.52, 1.0],
        "midnightblue" => [0.1, 0.1, 0.44, 1.0],
        "mintcream" => [0.96, 1.0, 0.98, 1.0],
        "mistyrose" => [1.0, 0.89, 0.88, 1.0],
        "moccasin" => [1.0, 0.89, 0.71, 1.0],
        "navajowhite" => [1.0, 0.87, 0.68, 1.0],
        "oldlace" => [0.99, 0.96, 0.9, 1.0],
        "olivedrab" => [0.42, 0.56, 0.14, 1.0],
        "orangered" => [1.0, 0.27, 0.0, 1.0],
        "orchid" => [0.85, 0.44, 0.84, 1.0],
        "palegoldenrod" => [0.93, 0.91, 0.67, 1.0],
        "palegreen" => [0.6, 0.98, 0.6, 1.0],
        "paleturquoise" => [0.69, 0.93, 0.93, 1.0],
        "palevioletred" => [0.86, 0.44, 0.58, 1.0],
        "papayawhip" => [1.0, 0.94, 0.84, 1.0],
        "peachpuff" => [1.0, 0.85, 0.73, 1.0],
        "peru" => [0.8, 0.52, 0.25, 1.0],
        "pink" => [1.0, 0.75, 0.8, 1.0],
        "plum" => [0.87, 0.63, 0.87, 1.0],
        "powderblue" => [0.69, 0.88, 0.9, 1.0],
        "rosybrown" => [0.74, 0.56, 0.56, 1.0],
        "royalblue" => [0.25, 0.41, 0.88, 1.0],
        "saddlebrown" => [0.55, 0.27, 0.07, 1.0],
        "salmon" => [0.98, 0.5, 0.45, 1.0],
        "sandybrown" => [0.96, 0.64, 0.38, 1.0],
        "seagreen" => [0.18, 0.55, 0.34, 1.0],
        "seashell" => [1.0, 0.96, 0.93, 1.0],
        "sienna" => [0.63, 0.32, 0.18, 1.0],
        "skyblue" => [0.53, 0.81, 0.92, 1.0],
        "slateblue" => [0.42, 0.35, 0.8, 1.0],
        "slategray" => [0.44, 0.5, 0.56, 1.0],
        "slategrey" => [0.44, 0.5, 0.56, 1.0],
        "snow" => [1.0, 0.98, 0.98, 1.0],
        "springgreen" => [0.0, 1.0, 0.5, 1.0],
        "steelblue" => [0.27, 0.51, 0.71, 1.0],
        "tan" => [0.82, 0.71, 0.55, 1.0],
        "thistle" => [0.85, 0.75, 0.85, 1.0],
        "tomato" => [1.0, 0.39, 0.28, 1.0],
        "turquoise" => [0.25, 0.88, 0.82, 1.0],
        "violet" => [0.93, 0.51, 0.93, 1.0],
        "wheat" => [0.96, 0.87, 0.7, 1.0],
        "whitesmoke" => [0.96, 0.96, 0.96, 1.0],
        "yellowgreen" => [0.6, 0.8, 0.2, 1.0],
        "transparent" => [0.0, 0.0, 0.0, 0.0]
    )
    
    return get(named_colors, lowercase(name), nothing)
end

# ////////////////////////////////////////////////////////////////////////
function parse_stroke_width(width_str::String)::Float64
    if isempty(width_str)
        return 1.0  # Default stroke width
    end
    width_str = replace(width_str, r"[a-z]+" => "")
    return parse(Float64, width_str)
end


# ////////////////////////////////////////////////////////////////////////
function apply_opacity(color::Vector{Float64}, opacity::Float64)::Vector{Float64}
    result = copy(color)
    result[4] *= opacity
    
    return result
end


# ////////////////////////////////////////////////////////////////////////
function parse_svg_file(filename::String; num_points::Int=24)::Dict{String, Any}
    doc = read(filename, XML.Node)
    svg_root = find_svg_root(doc)
    svg_tree = convert_svg_to_paths(svg_root, num_points=num_points)
    return clean_svg_tree(svg_tree)
end


# ////////////////////////////////////////////////////////////////////////
function find_svg_root(doc::XML.Node)::XML.Node
    for child in XML.children(doc)
        if XML.nodetype(child) == XML.Element && XML.tag(child) == "svg"
            return child
        end
    end
    for child in XML.children(doc)
        if XML.nodetype(child) == XML.Element
            for grandchild in XML.children(child)
                if XML.nodetype(grandchild) == XML.Element && XML.tag(grandchild) == "svg"
                    return grandchild
                end
            end
        end
    end
    
    error("No SVG element found in the document")
end


# ////////////////////////////////////////////////////////////////////////
function convert_svg_to_paths(svg_node::XML.Node; num_points::Int=24, parent_transform::Matrix{Float64}=[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])::Union{Dict{String, Any}, Nothing}
    attrs = get_attributes(svg_node)
    id = get(attrs, "id", "")
    tag = XML.tag(svg_node)
    
    if tag == "g" || tag == "svg"
        children = Dict{String, Any}[]
        element_transform = get_transform_matrix(attrs)
        combined_transform = parent_transform * element_transform
        for child in XML.children(svg_node)
            if XML.nodetype(child) == XML.Element
                child_node = convert_svg_to_paths(child, num_points=num_points, parent_transform=combined_transform)
                if child_node !== nothing
                    push!(children, child_node)
                end
            end
        end
        return Dict("type" => "group", "id" => id, "children" => children)
    elseif tag in ["rect", "circle", "ellipse", "line", "polyline", "polygon", "path"]
        element_transform = get_transform_matrix(attrs)
        combined_transform = parent_transform * element_transform
        path = shape_to_path(svg_node, num_points)
        transformed_subpaths = [apply_transform_to_subpath(sp, combined_transform) for sp in path["subpaths"]]
        
        return Dict("type" => "path", "id" => id, "style" => path["style"], "subpaths" => transformed_subpaths)
    else
        return nothing
    end
end


# ////////////////////////////////////////////////////////////////////////
function apply_transform_to_subpath(subpath::Dict{String, Any}, transform_matrix::Matrix{Float64})::Dict{String, Any}
    if isapprox(transform_matrix, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
        return subpath
    end
    points = subpath["points"]
    transformed_points = Vector{Vector{Float64}}(undef, length(points))
    for i in 1:length(points)
        p = points[i]
        x = transform_matrix[1,1] * p[1] + transform_matrix[1,2] * p[2] + transform_matrix[1,3]
        y = transform_matrix[2,1] * p[1] + transform_matrix[2,2] * p[2] + transform_matrix[2,3]
        
        transformed_points[i] = [x, y]
    end
    return Dict("points" => transformed_points, "closed" => subpath["closed"])
end


# ////////////////////////////////////////////////////////////////////////
function clean_svg_tree(node::Union{Dict{String, Any}, Nothing})::Union{Dict{String, Any}, Nothing}
    if node === nothing
        return nothing
    end
    
    if node["type"] == "group"
        cleaned_children = Dict{String, Any}[]
        
        for child in node["children"]
            cleaned_child = clean_svg_tree(child)
            if cleaned_child !== nothing
                push!(cleaned_children, cleaned_child)
            end
        end
        return Dict("type" => "group", "id" => node["id"], "children" => cleaned_children)
    elseif node["type"] == "path"
        return node
    else
        return nothing
    end
end


# ////////////////////////////////////////////////////////////////////////
function get_attributes(node::XML.Node)::Dict{String, String}
    attrs = Dict{String, String}()
    
    if XML.nodetype(node) == XML.Element && !isnothing(XML.attributes(node))
        for (key, value) in XML.attributes(node)
            attrs[key] = value
        end
    end
    if haskey(attrs, "style")
        style_attrs = parse_style_attribute(attrs["style"])
        for (key, value) in style_attrs
            attrs[key] = value
        end
    end
    
    return attrs
end


# ////////////////////////////////////////////////////////////////////////
function parse_style_attribute(style_str::String)::Dict{String, String}
    style_attrs = Dict{String, String}()
    for style_part in Base.split(style_str, ";")
        if isempty(style_part)
            continue
        end
        key_value = Base.split(style_part, ":")
        if length(key_value) == 2
            key = strip(key_value[1])
            value = strip(key_value[2])
            style_attrs[key] = value
        end
    end
    
    return style_attrs
end


# ////////////////////////////////////////////////////////////////////////
function get_transform_matrix(attrs::Dict{String, String})::Matrix{Float64}
    matrix = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    if haskey(attrs, "transform")
        transform_str = attrs["transform"]
        matrix = parse_transform(transform_str)
    end
    
    return matrix
end


# ////////////////////////////////////////////////////////////////////////
function parse_transform(transform_str::String)::Matrix{Float64}
    matrix = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    transform_functions = extract_transform_functions(transform_str)
    for transform_func in reverse(transform_functions)
        matrix = apply_transform_function(transform_func, matrix)
    end
    
    return matrix
end


# ////////////////////////////////////////////////////////////////////////
function apply_transform_function(transform_func::String, matrix::Matrix{Float64})
    translate_regex = r"translate\(\s*([^,\)]+)(?:,\s*([^,\)]+))?\s*\)"
    scale_regex = r"scale\(\s*([^,\)]+)(?:,\s*([^,\)]+))?\s*\)"
    rotate_regex = r"rotate\(\s*([^,\)]+)(?:,\s*([^,\)]+)(?:,\s*([^,\)]+))?)?\s*\)"
    matrix_regex = r"matrix\(\s*([^,\s]+)[\s,]+([^,\s]+)[\s,]+([^,\s]+)[\s,]+([^,\s]+)[\s,]+([^,\s]+)[\s,]+([^,\s]+)\s*\)"
    skewX_regex = r"skewX\(\s*([^,\)]+)\s*\)"
    skewY_regex = r"skewY\(\s*([^,\)]+)\s*\)"
    m = match(translate_regex, transform_func)
    if !isnothing(m)
        tx = parse(Float64, m.captures[1])
        ty = isnothing(m.captures[2]) ? 0.0 : parse(Float64, m.captures[2])
        translate_matrix = [1.0 0.0 tx; 0.0 1.0 ty; 0.0 0.0 1.0]
        return translate_matrix * matrix
    end
    m = match(scale_regex, transform_func)
    if !isnothing(m)
        sx = parse(Float64, m.captures[1])
        sy = isnothing(m.captures[2]) ? sx : parse(Float64, m.captures[2])
        scale_matrix = [sx 0.0 0.0; 0.0 sy 0.0; 0.0 0.0 1.0]
        return scale_matrix * matrix
    end
    m = match(rotate_regex, transform_func)
    if !isnothing(m)
        angle_deg = parse(Float64, m.captures[1])
        angle_rad = deg2rad(angle_deg)
        cos_a = cos(angle_rad)
        sin_a = sin(angle_rad)
        
        if isnothing(m.captures[2]) || isnothing(m.captures[3])
            rotate_matrix = [cos_a -sin_a 0.0; sin_a cos_a 0.0; 0.0 0.0 1.0]
            return rotate_matrix * matrix
        else
            cx = parse(Float64, m.captures[2])
            cy = parse(Float64, m.captures[3])
            translate_to_origin = [1.0 0.0 -cx; 0.0 1.0 -cy; 0.0 0.0 1.0]
            rotate_matrix = [cos_a -sin_a 0.0; sin_a cos_a 0.0; 0.0 0.0 1.0]
            translate_back = [1.0 0.0 cx; 0.0 1.0 cy; 0.0 0.0 1.0]
            
            transform_matrix = translate_back * rotate_matrix * translate_to_origin
            return transform_matrix * matrix
        end
    end
    m = match(matrix_regex, transform_func)
    if !isnothing(m)
        a = parse(Float64, m.captures[1])
        b = parse(Float64, m.captures[2])
        c = parse(Float64, m.captures[3])
        d = parse(Float64, m.captures[4])
        e = parse(Float64, m.captures[5])
        f = parse(Float64, m.captures[6])
        
        matrix_transform = [a c e; b d f; 0.0 0.0 1.0]
        return matrix_transform * matrix
    end
    m = match(skewX_regex, transform_func)
    if !isnothing(m)
        angle_deg = parse(Float64, m.captures[1])
        angle_rad = deg2rad(angle_deg)
        tan_a = tan(angle_rad)
        
        skew_matrix = [1.0 tan_a 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        return skew_matrix * matrix
    end
    m = match(skewY_regex, transform_func)
    if !isnothing(m)
        angle_deg = parse(Float64, m.captures[1])
        angle_rad = deg2rad(angle_deg)
        tan_a = tan(angle_rad)
        
        skew_matrix = [1.0 0.0 0.0; tan_a 1.0 0.0; 0.0 0.0 1.0]
        return skew_matrix * matrix
    end
    return matrix
end


# ////////////////////////////////////////////////////////////////////////
function extract_transform_functions(transform_str::String)
    transform_functions = String[]
    transform_regex = r"([a-zA-Z]+\([^)]*\))"
    for m in eachmatch(transform_regex, transform_str)
        push!(transform_functions, m.captures[1])
    end
    
    return transform_functions
end


# ////////////////////////////////////////////////////////////////////////
function matrix_to_transform_str(matrix::Matrix{Float64})
    a = matrix[1, 1]
    b = matrix[2, 1]
    c = matrix[1, 2]
    d = matrix[2, 2]
    e = matrix[1, 3]
    f = matrix[2, 3]
    return @sprintf("matrix(%.6f,%.6f,%.6f,%.6f,%.6f,%.6f)", a, b, c, d, e, f)
end


# ////////////////////////////////////////////////////////////////////////
function shape_to_path(svg_node::XML.Node, num_points::Int)
    tag = XML.tag(svg_node)
    attrs = get_attributes(svg_node)
    style = extract_style_attributes(attrs)
    if tag == "rect"
        return rect_to_path(attrs, style)
    elseif tag == "circle"
        return circle_to_path(attrs, style, num_points)
    elseif tag == "ellipse"
        return ellipse_to_path(attrs, style, num_points)
    elseif tag == "line"
        return line_to_path(attrs, style)
    elseif tag == "polyline"
        return parse_path_points(attrs, style, false)
    elseif tag == "polygon"
        return parse_path_points(attrs, style, true)
    elseif tag == "path"
        return path_element_to_path(attrs, style, num_points)
    else
        error("Unsupported shape type: $tag")
    end
end


# ////////////////////////////////////////////////////////////////////////
function extract_style_attributes(attrs::Dict{String, String})
    fill_color = get(attrs, "fill", "black")
    stroke_color = get(attrs, "stroke", "none")
    stroke_width_str = get(attrs, "stroke-width", "1")
    opacity_str = get(attrs, "opacity", "1")
    fill_opacity_str = get(attrs, "fill-opacity", "1")
    stroke_opacity_str = get(attrs, "stroke-opacity", "1")
    fill_rgba = parse_color(fill_color)
    stroke_rgba = parse_color(stroke_color)
    opacity = parse(Float64, opacity_str)
    fill_opacity = parse(Float64, fill_opacity_str)
    stroke_opacity = parse(Float64, stroke_opacity_str)
    fill_rgba[4] *= fill_opacity * opacity
    stroke_rgba[4] *= stroke_opacity * opacity
    stroke_width = parse_stroke_width(stroke_width_str)
    other_style = Dict{String, String}()
    for (key, value) in attrs
        if key != "fill" && key != "stroke" && key != "stroke-width" && 
           key != "opacity" && key != "fill-opacity" && key != "stroke-opacity"
            other_style[key] = value
        end
    end
    style = Dict{String, Any}(
        "fill" => fill_rgba,
        "stroke" => stroke_rgba,
        "stroke_width" => stroke_width,
        "opacity" => opacity
    )
    for (key, value) in other_style
        style[key] = value
    end
    
    return style
end


# ////////////////////////////////////////////////////////////////////////
function rect_to_path(attrs::Dict{String, String}, style::Dict{String, Any})::Dict{String, Any}
    x = parse(Float64, get(attrs, "x", "0"))
    y = parse(Float64, get(attrs, "y", "0"))
    width = parse(Float64, attrs["width"])
    height = parse(Float64, attrs["height"])
    points = [
        [x, y],
        [x + width, y],
        [x + width, y + height],
        [x, y + height]
    ]
    subpath = Dict("points" => points, "closed" => true)
    return Dict("style" => style, "subpaths" => [subpath])
end


# ////////////////////////////////////////////////////////////////////////
function circle_to_path(attrs::Dict{String, String}, style::Dict{String, Any}, num_points::Int)::Dict{String, Any}
    cx = parse(Float64, get(attrs, "cx", "0"))
    cy = parse(Float64, get(attrs, "cy", "0"))
    r = parse(Float64, attrs["r"])
    points = Vector{Vector{Float64}}(undef, num_points)
    for i in 1:num_points
        angle = 2π * (i - 1) / num_points
        x = cx + r * cos(angle)
        y = cy + r * sin(angle)
        points[i] = [x, y]
    end
    subpath = Dict("points" => points, "closed" => true)
    return Dict("style" => style, "subpaths" => [subpath])
end


# ////////////////////////////////////////////////////////////////////////
function ellipse_to_path(attrs::Dict{String, String}, style::Dict{String, Any}, num_points::Int)::Dict{String, Any}
    cx = parse(Float64, get(attrs, "cx", "0"))
    cy = parse(Float64, get(attrs, "cy", "0"))
    rx = parse(Float64, attrs["rx"])
    ry = parse(Float64, attrs["ry"])
    points = Vector{Vector{Float64}}(undef, num_points)
    for i in 1:num_points
        angle = 2π * (i - 1) / num_points
        x = cx + rx * cos(angle)
        y = cy + ry * sin(angle)
        points[i] = [x, y]
    end
    subpath = Dict("points" => points, "closed" => true)
    return Dict("style" => style, "subpaths" => [subpath])
end


# ////////////////////////////////////////////////////////////////////////
function line_to_path(attrs::Dict{String, String}, style::Dict{String, Any})::Dict{String, Any}
    x1 = parse(Float64, get(attrs, "x1", "0"))
    y1 = parse(Float64, get(attrs, "y1", "0"))
    x2 = parse(Float64, get(attrs, "x2", "0"))
    y2 = parse(Float64, get(attrs, "y2", "0"))
    points = [[x1, y1], [x2, y2]]
    subpath = Dict("points" => points, "closed" => false)
    return Dict("style" => style, "subpaths" => [subpath])
end


# ////////////////////////////////////////////////////////////////////////
function parse_path_points(attrs::Dict{String, String}, style::Dict{String, Any}, closed::Bool)::Dict{String, Any}
    points_str = attrs["points"]
    points = Vector{Float64}[]
    for point_str in Base.split(points_str)
        if isempty(point_str)
            continue
        end
        
        coords = Base.split(point_str, ',')
        if length(coords) >= 2
            x = parse(Float64, coords[1])
            y = parse(Float64, coords[2])
            push!(points, [x, y])
        end
    end
    subpath = Dict("points" => points, "closed" => closed)
    return Dict("style" => style, "subpaths" => [subpath])
end


# ////////////////////////////////////////////////////////////////////////
function path_element_to_path(attrs::Dict{String, String}, style::Dict{String, Any}, num_points::Int)::Union{Dict{String, Any}, Vector{Dict{String, Any}}}
    path_data = attrs["d"]
    return parse_path_data(path_data, style, num_points)
end


# ////////////////////////////////////////////////////////////////////////
function parse_path_data(path_data::String, style::Dict{String, Any}, num_points::Int)::Union{Dict{String, Any}, Vector{Dict{String, Any}}}
    tokens = tokenize_path_data(path_data)
    subpaths = Dict{String, Any}[]
    points = Vector{Float64}[]
    current_point = [0.0, 0.0]
    start_point = [0.0, 0.0]
    closed = false
    first_move = true
    
    i = 1
    while i <= length(tokens)
        command = tokens[i]
        i += 1
        
        if command == "M" || command == "m"
            if !first_move && !isempty(points)
                push!(subpaths, Dict("points" => points, "closed" => closed))
                points = Vector{Float64}[]
                closed = false
            end
            first_move = false
            relative = command == "m"
            x = parse(Float64, tokens[i])
            y = parse(Float64, tokens[i+1])
            i += 2
            
            if relative
                current_point = [current_point[1] + x, current_point[2] + y]
            else
                current_point = [x, y]
            end
            
            start_point = copy(current_point)
            push!(points, copy(current_point))
            
        elseif command == "L" || command == "l"
            relative = command == "l"
            x = parse(Float64, tokens[i])
            y = parse(Float64, tokens[i+1])
            i += 2
            
            if relative
                current_point = [current_point[1] + x, current_point[2] + y]
            else
                current_point = [x, y]
            end
            
            push!(points, copy(current_point))
            
        elseif command == "H" || command == "h"
            relative = command == "h"
            x = parse(Float64, tokens[i])
            i += 1
            
            if relative
                current_point = [current_point[1] + x, current_point[2]]
            else
                current_point = [x, current_point[2]]
            end
            
            push!(points, copy(current_point))
            
        elseif command == "V" || command == "v"
            relative = command == "v"
            y = parse(Float64, tokens[i])
            i += 1
            
            if relative
                current_point = [current_point[1], current_point[2] + y]
            else
                current_point = [current_point[1], y]
            end
            
            push!(points, copy(current_point))
            
        elseif command == "C" || command == "c"
            relative = command == "c"
            
            x1 = parse(Float64, tokens[i])
            y1 = parse(Float64, tokens[i+1])
            x2 = parse(Float64, tokens[i+2])
            y2 = parse(Float64, tokens[i+3])
            x = parse(Float64, tokens[i+4])
            y = parse(Float64, tokens[i+5])
            i += 6
            
            if relative
                x1 += current_point[1]
                y1 += current_point[2]
                x2 += current_point[1]
                y2 += current_point[2]
                x += current_point[1]
                y += current_point[2]
            end
            for t in range(0, 1, length=num_points)
                bx = (1-t)^3 * current_point[1] + 3*(1-t)^2*t * x1 + 3*(1-t)*t^2 * x2 + t^3 * x
                by = (1-t)^3 * current_point[2] + 3*(1-t)^2*t * y1 + 3*(1-t)*t^2 * y2 + t^3 * y
                
                push!(points, [bx, by])
            end
            
            current_point = [x, y]
            
        elseif command == "S" || command == "s"
            relative = command == "s"
            if length(points) >= 2
                prev_control = points[end-1]
                x1 = 2 * current_point[1] - prev_control[1]
                y1 = 2 * current_point[2] - prev_control[2]
            else
                x1 = current_point[1]
                y1 = current_point[2]
            end
            
            x2 = parse(Float64, tokens[i])
            y2 = parse(Float64, tokens[i+1])
            x = parse(Float64, tokens[i+2])
            y = parse(Float64, tokens[i+3])
            i += 4
            
            if relative
                x2 += current_point[1]
                y2 += current_point[2]
                x += current_point[1]
                y += current_point[2]
            end
            for t in range(0, 1, length=num_points)
                bx = (1-t)^3 * current_point[1] + 3*(1-t)^2*t * x1 + 3*(1-t)*t^2 * x2 + t^3 * x
                by = (1-t)^3 * current_point[2] + 3*(1-t)^2*t * y1 + 3*(1-t)*t^2 * y2 + t^3 * y
                
                push!(points, [bx, by])
            end
            
            current_point = [x, y]
            
        elseif command == "Q" || command == "q"
            relative = command == "q"
            
            x1 = parse(Float64, tokens[i])
            y1 = parse(Float64, tokens[i+1])
            x = parse(Float64, tokens[i+2])
            y = parse(Float64, tokens[i+3])
            i += 4
            
            if relative
                x1 += current_point[1]
                y1 += current_point[2]
                x += current_point[1]
                y += current_point[2]
            end
            for t in range(0, 1, length=num_points)
                bx = (1-t)^2 * current_point[1] + 2*(1-t)*t * x1 + t^2 * x
                by = (1-t)^2 * current_point[2] + 2*(1-t)*t * y1 + t^2 * y
                
                push!(points, [bx, by])
            end
            
            current_point = [x, y]
            
        elseif command == "T" || command == "t"
            relative = command == "t"
            if length(points) >= 2
                prev_control = points[end-1]
                x1 = 2 * current_point[1] - prev_control[1]
                y1 = 2 * current_point[2] - prev_control[2]
            else
                x1 = current_point[1]
                y1 = current_point[2]
            end
            
            x = parse(Float64, tokens[i])
            y = parse(Float64, tokens[i+1])
            i += 2
            
            if relative
                x += current_point[1]
                y += current_point[2]
            end
            for t in range(0, 1, length=num_points)
                bx = (1-t)^2 * current_point[1] + 2*(1-t)*t * x1 + t^2 * x
                by = (1-t)^2 * current_point[2] + 2*(1-t)*t * y1 + t^2 * y
                
                push!(points, [bx, by])
            end
            
            current_point = [x, y]
            
        elseif command == "A" || command == "a"
            relative = command == "a"
            
            rx = parse(Float64, tokens[i])
            ry = parse(Float64, tokens[i+1])
            x_axis_rotation = parse(Float64, tokens[i+2])
            large_arc_flag = parse(Float64, tokens[i+3])
            sweep_flag = parse(Float64, tokens[i+4])
            x = parse(Float64, tokens[i+5])
            y = parse(Float64, tokens[i+6])
            i += 7
            
            if relative
                x += current_point[1]
                y += current_point[2]
            end
            arc_points = elliptical_arc_to_points(
                current_point[1], current_point[2],
                rx, ry,
                deg2rad(x_axis_rotation),
                large_arc_flag != 0,
                sweep_flag != 0,
                x, y,
                num_points
            )
            
            append!(points, arc_points)
            current_point = [x, y]
            
        elseif command == "Z" || command == "z"
            if length(points) > 0
                if points[end] != start_point
                    push!(points, copy(start_point))
                end
                closed = true
            end
            
        else
            @warn "Unknown path command: $command"
        end
    end
    if !isempty(points)
        push!(subpaths, Dict("points" => points, "closed" => closed))
    end
    return Dict("style" => style, "subpaths" => subpaths)
end


# ////////////////////////////////////////////////////////////////////////
function tokenize_path_data(path_data::String)::Vector{SubString{String}}
    path_data = replace(path_data, "," => " ")
    for cmd in "MmLlHhVvCcSsQqTtAaZz"
        path_data = replace(path_data, cmd => " $cmd ")
    end
    tokens = filter(!isempty, Base.split(path_data))
    
    return tokens
end


# ////////////////////////////////////////////////////////////////////////
function elliptical_arc_to_points(x1, y1, rx, ry, phi, large_arc, sweep, x2, y2, num_points)::Vector{Vector{Float64}}
    rx = abs(rx)
    ry = abs(ry)
    if rx < 1e-10 || ry < 1e-10
        return [[x2, y2]]
    end
    cos_phi = cos(phi)
    sin_phi = sin(phi)
    dx = (x1 - x2) / 2
    dy = (y1 - y2) / 2
    x1p = cos_phi * dx + sin_phi * dy
    y1p = -sin_phi * dx + cos_phi * dy
    lambda = (x1p / rx)^2 + (y1p / ry)^2
    if lambda > 1
        rx *= sqrt(lambda)
        ry *= sqrt(lambda)
    end
    sign = large_arc == sweep ? -1 : 1
    sq = ((rx*ry)^2 - (rx*y1p)^2 - (ry*x1p)^2) / ((rx*y1p)^2 + (ry*x1p)^2)
    sq = max(0, sq)  # Ensure non-negative
    coef = sign * sqrt(sq)
    
    cxp = coef * (rx * y1p) / ry
    cyp = coef * (-ry * x1p) / rx
    cx = cos_phi * cxp - sin_phi * cyp + (x1 + x2) / 2
    cy = sin_phi * cxp + cos_phi * cyp + (y1 + y2) / 2
    ux = (x1p - cxp) / rx
    uy = (y1p - cyp) / ry
    vx = (-x1p - cxp) / rx
    vy = (-y1p - cyp) / ry
    theta1 = angle_between(1.0, 0.0, ux, uy)
    delta = angle_between(ux, uy, vx, vy)
    if !sweep && delta > 0
        delta -= 2π
    elseif sweep && delta < 0
        delta += 2π
    end
    points = Vector{Float64}[]
    for i in 0:num_points
        t = theta1 + delta * i / num_points
        px = cos(t)
        py = sin(t)
        x = cos_phi * rx * px - sin_phi * ry * py + cx
        y = sin_phi * rx * px + cos_phi * ry * py + cy
        
        push!(points, [x, y])
    end
    
    return points
end


# ////////////////////////////////////////////////////////////////////////
function angle_between(ux, uy, vx, vy)::Float64
    dot = ux * vx + uy * vy
    len_u = sqrt(ux^2 + uy^2)
    len_v = sqrt(vx^2 + vy^2)
    angle = acos(clamp(dot / (len_u * len_v), -1.0, 1.0))
    if ux * vy - uy * vx < 0
        angle = -angle
    end
    
    return angle
end


# ////////////////////////////////////////////////////////////////////////
function collect_paths(node::Dict{String, Any})::Vector{Dict{String, Any}}
    paths = Dict{String, Any}[]
    if node["type"] == "group"
        for child in node["children"]
            append!(paths, collect_paths(child))
        end
    elseif node["type"] == "path"
        push!(paths, Dict("style" => node["style"], "subpaths" => node["subpaths"]))
    end
    return paths
end
export collect_paths


# ////////////////////////////////////////////////////////////////////////
function paths_to_edges(paths::Vector{Dict{String, Any}})::Tuple{Vector{Tuple{Float64, Float64}}, Vector{Vector{Vector{Tuple{Int, Int}}}}}
    all_vertices = Vector{Tuple{Float64, Float64}}()
    all_paths_edges = Vector{Vector{Vector{Tuple{Int, Int}}}}()
    vertex_map = Dict{Tuple{Float64, Float64}, Int}()
    for path in paths
        path_subpaths_edges = Vector{Vector{Tuple{Int, Int}}}()
        for subpath in path["subpaths"]
            subpath_edges = Vector{Tuple{Int, Int}}()
            points = subpath["points"]
            for i in 1:length(points)-1
                p1 = (points[i][1], points[i][2])
                p2 = (points[i+1][1], points[i+1][2])
                if !haskey(vertex_map, p1)
                    push!(all_vertices, p1)
                    vertex_map[p1] = length(all_vertices)
                end
                if !haskey(vertex_map, p2)
                    push!(all_vertices, p2)
                    vertex_map[p2] = length(all_vertices)
                end
                push!(subpath_edges, (vertex_map[p1], vertex_map[p2]))
            end
            push!(path_subpaths_edges, subpath_edges)
        end
        push!(all_paths_edges, path_subpaths_edges)
    end
    return (all_vertices, all_paths_edges)
end
export paths_to_edges




# ////////////////////////////////////////////////////////////////////////
function svg_to_hpc_edges(node)
    paths = collect_paths(node)
    v = []
    for path in paths
        style = path["style"]
        line_color = Point4d(style["stroke"]...)
        line_width = max(1, Int(round(style["stroke_width"])))
        face_color = Point4d(style["fill"]...)
        for subpath in path["subpaths"]
            points = [[p[1], -p[2]] for p in subpath["points"]]
            if subpath["closed"]
                if length(points) > 1 && points[1] != points[end]
                    push!(points, points[1])
                end
            end
            geometry = Geometry()
            for i in 1:length(points)-1
                if points[i] != points[i+1]
                    addHull(geometry, [points[i], points[i+1]])
                end
            end
            result = PROPERTIES(Hpc(geometry), Properties("line_color" => line_color, "line_width" => line_width, "face_color" => face_color))
            push!(v, result)
        end
    end
    return length(v) > 0 ? STRUCT(v) : Hpc()
end

# ////////////////////////////////////////////////////////////////////////
function svg_to_lar_edges(node)
    paths = collect_paths(node)
    
    # Collect all vertices and edges
    all_vertices = Vector{Vector{Float64}}()
    all_edges = Vector{Tuple{Int, Int}}()
    vertex_map = Dict{Tuple{Float64, Float64}, Int}()
    
    for path in paths
        for subpath in path["subpaths"]
            points = [[p[1], -p[2]] for p in subpath["points"]]
            if subpath["closed"]
                if length(points) > 1 && points[1] != points[end]
                    push!(points, points[1])
                end
            end
            
            # Add vertices and edges
            for i in 1:length(points)-1
                if points[i] != points[i+1]
                    p1 = (points[i][1], points[i][2])
                    p2 = (points[i+1][1], points[i+1][2])
                    
                    # Get or create vertex index for p1
                    if !haskey(vertex_map, p1)
                        push!(all_vertices, [p1[1], p1[2]])
                        vertex_map[p1] = length(all_vertices)
                    end
                    
                    # Get or create vertex index for p2
                    if !haskey(vertex_map, p2)
                        push!(all_vertices, [p2[1], p2[2]])
                        vertex_map[p2] = length(all_vertices)
                    end
                    
                    # Add edge
                    push!(all_edges, (vertex_map[p1], vertex_map[p2]))
                end
            end
        end
    end
    
    # Convert vertices to matrix format (points as columns)
    if isempty(all_vertices)
        V = Matrix{Float64}(undef, 2, 0)
    else
        V = hcat(all_vertices...)
    end
    
    # Convert edges to Cells format
    EV = [[e[1], e[2]] for e in all_edges]
    
    # Create and return Lar object
    return Lar(V, Dict{Symbol, Cells}(:EV => EV))
end


# ////////////////////////////////////////////////////////////////////////////
function build_output(V::Matrix{Float64}, FV::Vector, output_type::Symbol)
    if output_type == :lar
        all_edges = Set{Tuple{Int, Int}}()
        for fv in FV
            push!(all_edges, (min(fv[1], fv[2]), max(fv[1], fv[2])))
            push!(all_edges, (min(fv[2], fv[3]), max(fv[2], fv[3])))
            push!(all_edges, (min(fv[3], fv[1]), max(fv[3], fv[1])))
        end
        EV = [[e[1], e[2]] for e in all_edges]
        return Lar(V, Dict{Symbol, Cells}(:EV => EV, :FV => FV))
    else
        geometry = Geometry()
        for fv in FV
            v1, v2, v3 = V[:, fv[1]], V[:, fv[2]], V[:, fv[3]]
            addHull(geometry, [v1, v2, v3])
        end
        return Hpc(geometry)
    end
end

# ////////////////////////////////////////////////////////////////////////////
export svg_to_triangles
function svg_to_triangles(svg; triangulate::Symbol=:global, output_type::Symbol=:hpc)
    
    println("svg_to_triangles triangulate=$(triangulate) output_type=$(output_type)")
    
    # Extract all paths and edges
    paths = collect_paths(svg)
    all_vertices, all_shapes = paths_to_edges(paths)
    
    # Flatten all shapes to get all edges
    all_edges_flat = vcat([vcat(shape...) for shape in all_shapes]...)
    V = hcat([[v[1], v[2]] for v in all_vertices]...)
    EV_matrix = hcat([[e[1], e[2]] for e in all_edges_flat]...)
    
    println("\nTotal vertices: $(size(V, 2))")
    println("Total shapes: $(length(all_shapes))")
    println("Total edges (with duplicates): $(length(all_edges_flat))")
    
    # Process based on classification mode
    @time begin
        if triangulate == :local
            # Triangulate and classify per-shape
            v = []
            
            for shape in all_shapes
                shape_edges_flat = vcat(shape...)
                isempty(shape_edges_flat) && continue
                
                # Get unique vertices and create local matrix
                shape_vertex_set = Set{Int}()
                for edge in shape_edges_flat
                    push!(shape_vertex_set, edge[1])
                    push!(shape_vertex_set, edge[2])
                end
                shape_vertices_list = sort(collect(shape_vertex_set))
                V_shape = hcat([[all_vertices[idx][1], all_vertices[idx][2]] for idx in shape_vertices_list]...)
                
                # Remap edges to local indices
                global_to_local = Dict(global_idx => local_idx for (local_idx, global_idx) in enumerate(shape_vertices_list))
                EV_shape = [[global_to_local[e[1]], global_to_local[e[2]]] for e in shape_edges_flat]
                
                # Triangulate and classify
                triin = Triangulate.TriangulateIO()
                triin.pointlist = V_shape
                triin.segmentlist = hcat(EV_shape...)
                (triout, _) = Triangulate.triangulate("pQ", triin)
                
                V_tri = triout.pointlist
                FV = [Int64.(triout.trianglelist[:, i]) for i in 1:size(triout.trianglelist, 2)]
                in_out = classify_triangles_robust(V_tri, FV, Matrix(V_shape'), EV_shape)
                FV_accepted = FV[findall(x -> x == 1, in_out)]
                
                if output_type == :lar
                    for fv in FV_accepted
                        push!(v, [shape_vertices_list[vx] for vx in fv])
                    end
                else
                    push!(v, build_output(V_tri, FV_accepted, :hpc))
                end
            end
            
            if output_type == :lar
                return build_output(V, v, :lar)
            else
                return STRUCT(v)
            end
        else
            # Triangulate globally (all shapes together)
            triin = Triangulate.TriangulateIO()
            triin.pointlist = V
            triin.segmentlist = EV_matrix
            (triout, _) = Triangulate.triangulate("pQ", triin)
            
            V_tri = triout.pointlist
            FV = [Int64.(triout.trianglelist[:, i]) for i in 1:size(triout.trianglelist, 2)]
            println("\nTriangles generated: $(length(FV))")
            
            # Classify triangles
            println("\nClassifying triangles...")
            V_orig_row = Matrix(V')
            EV_cells = [[e[1], e[2]] for e in all_edges_flat]
            in_out = classify_triangles_robust(V_tri, FV, V_orig_row, EV_cells)
            FV_accepted = FV[findall(x -> x == 1, in_out)]
            println("\nAccepted triangles: $(length(FV_accepted))")
            
            return build_output(V_tri, FV_accepted, output_type)
        end
    end
end

# ////////////////////////////////////////////////////////////////////////////
export svg_to_plasm
function svg_to_plasm(node; triangulate::Union{Nothing, Symbol}=nothing, output_type::Symbol=:hpc)
    if isnothing(triangulate)
        if output_type == :hpc
            return svg_to_hpc_edges(node)
        else
            return svg_to_lar_edges(node)
        end
    else
        return svg_to_triangles(node, triangulate=triangulate, output_type=output_type)
    end
end

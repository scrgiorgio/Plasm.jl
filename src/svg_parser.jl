

using XML
using LinearAlgebra
using Printf

export parse_svg_file, convert_svg_to_polylines, SVGPoint, SVGPoints, SVGNode, SVGGroup, SVGPolyline, Style

"""
Color parsing module for SVG files

This module provides functions to parse SVG color strings (like "rgb(255,0,0)" or "#FF0000")
to RGBA arrays [r, g, b, a] with values in the range 0.0-1.0.
"""

"""
    parse_color(color_str::String)

Parse an SVG color string to an RGBA array [r, g, b, a] with values in the range 0.0-1.0.
Supports hex colors, rgb(), rgba(), named colors, and "none".
"""
function parse_color(color_str::String)
    # Default color (black, fully opaque)
    rgba = [0.0, 0.0, 0.0, 1.0]
    
    # Handle "none" or empty string
    if isempty(color_str) || lowercase(color_str) == "none"
        return [0.0, 0.0, 0.0, 0.0]  # Transparent
    end
    
    # Handle hex colors
    if startswith(color_str, "#")
        return parse_hex_color(color_str)
    end
    
    # Handle rgb() and rgba() formats
    if startswith(lowercase(color_str), "rgb")
        return parse_rgb_color(color_str)
    end
    
    # Handle named colors
    named_color = get_named_color(color_str)
    if !isnothing(named_color)
        return named_color
    end
    
    # Return default if color format is not recognized
    return rgba
end

"""
    parse_hex_color(hex_str::String)

Parse a hex color string (like "#FF0000" or "#F00") to an RGBA array.
"""
function parse_hex_color(hex_str::String)
    hex = replace(hex_str, "#" => "")
    
    # Handle shorthand hex (#RGB)
    if length(hex) == 3
        r = parse(Int, hex[1] * hex[1], base=16) / 255.0
        g = parse(Int, hex[2] * hex[2], base=16) / 255.0
        b = parse(Int, hex[3] * hex[3], base=16) / 255.0
        return [r, g, b, 1.0]
    end
    
    # Handle shorthand hex with alpha (#RGBA)
    if length(hex) == 4
        r = parse(Int, hex[1] * hex[1], base=16) / 255.0
        g = parse(Int, hex[2] * hex[2], base=16) / 255.0
        b = parse(Int, hex[3] * hex[3], base=16) / 255.0
        a = parse(Int, hex[4] * hex[4], base=16) / 255.0
        return [r, g, b, a]
    end
    
    # Handle full hex (#RRGGBB)
    if length(hex) == 6
        r = parse(Int, hex[1:2], base=16) / 255.0
        g = parse(Int, hex[3:4], base=16) / 255.0
        b = parse(Int, hex[5:6], base=16) / 255.0
        return [r, g, b, 1.0]
    end
    
    # Handle full hex with alpha (#RRGGBBAA)
    if length(hex) == 8
        r = parse(Int, hex[1:2], base=16) / 255.0
        g = parse(Int, hex[3:4], base=16) / 255.0
        b = parse(Int, hex[5:6], base=16) / 255.0
        a = parse(Int, hex[7:8], base=16) / 255.0
        return [r, g, b, a]
    end
    
    # Default to black if format is not recognized
    return [0.0, 0.0, 0.0, 1.0]
end

"""
    parse_rgb_color(rgb_str::String)

Parse an rgb() or rgba() color string to an RGBA array.
"""
function parse_rgb_color(rgb_str::String)
    # Extract values from rgb(r,g,b) or rgba(r,g,b,a)
    is_rgba = startswith(lowercase(rgb_str), "rgba")
    
    # Extract the values inside the parentheses
    values_str = match(r"\((.*)\)", rgb_str).captures[1]
    values = Base.split(values_str, r"[,\s]+")
    
    # Parse RGB values
    r = parse_color_component(String(strip(values[1])))
    g = parse_color_component(String(strip(values[2])))
    b = parse_color_component(String(strip(values[3])))
    
    # Parse alpha value if present
    a = 1.0
    if is_rgba && length(values) >= 4
        a = parse(Float64, strip(values[4]))
    end
    
    return [r, g, b, a]
end

"""
    parse_color_component(component_str::AbstractString)

Parse a color component string, handling percentages and values.
"""
function parse_color_component(component_str::AbstractString)
    component_str = strip(component_str)
    
    # Handle percentage values
    if endswith(component_str, "%")
        percentage = parse(Float64, replace(component_str, "%" => ""))
        return clamp(percentage / 100.0, 0.0, 1.0)
    end
    
    # Handle numeric values (0-255)
    value = parse(Float64, component_str)
    if value <= 1.0 && !contains(component_str, ".")
        # Assume 0-255 range if value is an integer <= 1.0
        return clamp(value / 255.0, 0.0, 1.0)
    else
        # Assume 0.0-1.0 range if value is a float
        return clamp(value, 0.0, 1.0)
    end
end

"""
    get_named_color(name::String)

Get the RGBA values for a named color.
"""
function get_named_color(name::String)
    # Common SVG named colors
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

"""
    parse_stroke_width(width_str::String)

Parse a stroke-width string to a float value.
"""
function parse_stroke_width(width_str::String)
    # Handle empty string
    if isempty(width_str)
        return 1.0  # Default stroke width
    end
    
    # Remove units (px, pt, etc.) and parse as float
    width_str = replace(width_str, r"[a-z]+" => "")
    return parse(Float64, width_str)
end

"""
    apply_opacity(color::Vector{Float64}, opacity::Float64)

Apply an opacity value to a color's alpha channel.
"""
function apply_opacity(color::Vector{Float64}, opacity::Float64)
    # Create a copy of the color
    result = copy(color)
    
    # Apply opacity to alpha channel
    result[4] *= opacity
    
    return result
end

# Test the color parsing functions
function test_color_parsing()
    # Test hex colors
    @assert parse_color("#FF0000") ≈ [1.0, 0.0, 0.0, 1.0]
    @assert parse_color("#F00") ≈ [1.0, 0.0, 0.0, 1.0]
    @assert parse_color("#FF0000FF") ≈ [1.0, 0.0, 0.0, 1.0]
    @assert parse_color("#F00F") ≈ [1.0, 0.0, 0.0, 1.0]
    
    # Test rgb() and rgba() colors
    @assert parse_color("rgb(255, 0, 0)") ≈ [1.0, 0.0, 0.0, 1.0]
    @assert parse_color("rgb(100%, 0%, 0%)") ≈ [1.0, 0.0, 0.0, 1.0]
    @assert parse_color("rgba(255, 0, 0, 0.5)") ≈ [1.0, 0.0, 0.0, 0.5]
    
    # Test named colors
    @assert parse_color("red") ≈ [1.0, 0.0, 0.0, 1.0]
    @assert parse_color("blue") ≈ [0.0, 0.0, 1.0, 1.0]
    @assert parse_color("transparent") ≈ [0.0, 0.0, 0.0, 0.0]
    
    # Test opacity application
    @assert apply_opacity([1.0, 0.0, 0.0, 1.0], 0.5) ≈ [1.0, 0.0, 0.0, 0.5]
    
    println("All color parsing tests passed!")
end


"""
    SVGPoint

A 2D point with x and y coordinates.
"""
struct SVGPoint
    x::Float64
    y::Float64
end

"""
    Style

A structure to hold style attributes with RGBA arrays for colors.
"""
struct Style
    fill::Vector{Float64}        # RGBA array [r, g, b, a] with values 0.0-1.0
    stroke::Vector{Float64}      # RGBA array [r, g, b, a] with values 0.0-1.0
    stroke_width::Float64        # Line width as float
    opacity::Float64             # Overall opacity as float
    other::Dict{String, String}  # Other style attributes
end

"""
    SVGPoints

A polyline represented as a sequence of points.
"""
struct SVGPoints
    points::Vector{SVGPoint}
    closed::Bool
    style::Style
end

"""
    SVGNode

A node in the SVG tree structure.
"""
abstract type SVGNode end

"""
    SVGGroup

A group node in the SVG tree structure.
"""
struct SVGGroup <: SVGNode
    id::String
    children::Vector{SVGNode}
end

"""
    SVGPolyline

A polyline/polygon node in the SVG tree structure.
"""
struct SVGPolyline <: SVGNode
    id::String
    polyline::SVGPoints
end

"""
    parse_svg_file(filename::String; num_points::Int=24)

Parse an SVG file and return the root node of the SVG tree.

# Arguments
- `filename::String`: Path to the SVG file
- `num_points::Int=24`: Number of points to sample for curved elements

# Returns
- `SVGNode`: The root node of the converted tree
"""
function parse_svg_file(filename::String; num_points::Int=24)
    # Parse the SVG file using XML.jl
    doc = read(filename, XML.Node)
    
    # Find the root SVG element
    svg_root = find_svg_root(doc)
    
    # Convert the SVG tree to our internal representation
    svg_tree = convert_svg_to_polylines(svg_root, num_points=num_points)
    
    # Clean up the tree structure
    return clean_svg_tree(svg_tree)
end

"""
    find_svg_root(doc::XML.Node)

Find the root SVG element in an XML document.
"""
function find_svg_root(doc::XML.Node)
    # Look for the SVG element in the document
    for child in XML.children(doc)
        if XML.nodetype(child) == XML.Element && XML.tag(child) == "svg"
            return child
        end
    end
    
    # If not found at the top level, search recursively
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

"""
    convert_svg_to_polylines(svg_node::XML.Node; num_points::Int=24, parent_transform::Matrix{Float64}=[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])

Convert an SVG element to our internal tree structure with polylines.

# Arguments
- `svg_node::XML.Node`: The SVG element to convert
- `num_points::Int=24`: Number of points to sample for curved elements
- `parent_transform::Matrix{Float64}`: Transformation matrix from parent elements

# Returns
- `SVGNode`: The root node of the converted tree
"""
function convert_svg_to_polylines(svg_node::XML.Node; num_points::Int=24, parent_transform::Matrix{Float64}=[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    # Get the ID of the element
    attrs = get_attributes(svg_node)
    id = get(attrs, "id", "")
    
    # Process based on element type
    tag = XML.tag(svg_node)
    
    if tag == "g" || tag == "svg"
        # Process group element
        children = SVGNode[]
        
        # Get transformation matrix for this group
        element_transform = get_transform_matrix(attrs)
        
        # Combine with parent transform
        # For SVG, transformations are applied in order from right to left
        # So the parent transform should be applied after the element transform
        combined_transform = parent_transform * element_transform
        
        # Process children
        for child in XML.children(svg_node)
            if XML.nodetype(child) == XML.Element
                # Convert child to our internal representation with combined transform
                child_node = convert_svg_to_polylines(child, num_points=num_points, parent_transform=combined_transform)
                
                # Add to children if not nothing
                if child_node !== nothing
                    push!(children, child_node)
                end
            end
        end
        
        # Return group node (without transformation)
        return SVGGroup(id, children)
    elseif tag in ["rect", "circle", "ellipse", "line", "polyline", "polygon", "path"]
        # Get the element's own transformation
        element_transform = get_transform_matrix(attrs)
        
        # Combine with parent transform
        # The parent transform should be applied after the element transform
        combined_transform = parent_transform * element_transform
        
        # Convert shape to polyline
        polyline = shape_to_polyline(svg_node, num_points)
        
        # Apply transformation to all points in the polyline
        transformed_polyline = apply_transform_to_polyline(polyline, combined_transform)
        
        # Return polyline node
        return SVGPolyline(id, transformed_polyline)
    else
        # Ignore other elements
        return nothing
    end
end

"""
    apply_transform_to_polyline(polyline::SVGPoints, transform_matrix::Matrix{Float64})

Apply a transformation matrix to all points in a polyline.
"""
function apply_transform_to_polyline(polyline::SVGPoints, transform_matrix::Matrix{Float64})
    # If identity matrix, return polyline as is
    if isapprox(transform_matrix, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
        return polyline
    end
    
    # Apply transformation to each point
    transformed_points = Vector{SVGPoint}(undef, length(polyline.points))
    for i in 1:length(polyline.points)
        p = polyline.points[i]
        
        # Apply transformation matrix to point
        x = transform_matrix[1,1] * p.x + transform_matrix[1,2] * p.y + transform_matrix[1,3]
        y = transform_matrix[2,1] * p.x + transform_matrix[2,2] * p.y + transform_matrix[2,3]
        
        transformed_points[i] = SVGPoint(x, y)
    end
    
    # Return transformed polyline
    return SVGPoints(transformed_points, polyline.closed, polyline.style)
end

"""
    clean_svg_tree(node::SVGNode)

Clean up the SVG tree structure to ensure:
1. Intermediate nodes are only grouping nodes without transformations
2. Leaf nodes are polylines/polygons with all transformations applied
"""
function clean_svg_tree(node::SVGNode)
    if node isa SVGGroup
        # Process children recursively
        cleaned_children = SVGNode[]
        
        for child in node.children
            cleaned_child = clean_svg_tree(child)
            if cleaned_child !== nothing
                push!(cleaned_children, cleaned_child)
            end
        end
        
        # Return cleaned group
        return SVGGroup(node.id, cleaned_children)
    elseif node isa SVGPolyline
        # Leaf nodes are already polylines with transformations applied
        return node
    else
        # Ignore other node types
        return nothing
    end
end

"""
    get_attributes(node::XML.Node)

Get the attributes of an XML node as a dictionary.
"""
function get_attributes(node::XML.Node)
    attrs = Dict{String, String}()
    
    if XML.nodetype(node) == XML.Element && !isnothing(XML.attributes(node))
        for (key, value) in XML.attributes(node)
            attrs[key] = value
        end
    end
    
    # Also check for style attribute and parse it
    if haskey(attrs, "style")
        style_attrs = parse_style_attribute(attrs["style"])
        for (key, value) in style_attrs
            attrs[key] = value
        end
    end
    
    return attrs
end

"""
    parse_style_attribute(style_str::String)

Parse an SVG style attribute string into a dictionary.
"""
function parse_style_attribute(style_str::String)
    style_attrs = Dict{String, String}()
    
    # Split by semicolons
    for style_part in Base.split(style_str, ";")
        if isempty(style_part)
            continue
        end
        
        # Split by colon
        key_value = Base.split(style_part, ":")
        if length(key_value) == 2
            key = strip(key_value[1])
            value = strip(key_value[2])
            style_attrs[key] = value
        end
    end
    
    return style_attrs
end

"""
    get_transform_matrix(attrs::Dict{String, String})

Get the transformation matrix from element attributes.
"""
function get_transform_matrix(attrs::Dict{String, String})
    # Start with identity matrix
    matrix = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    
    # Check if transform attribute exists
    if haskey(attrs, "transform")
        # Parse transform string
        transform_str = attrs["transform"]
        matrix = parse_transform(transform_str)
    end
    
    return matrix
end

"""
    parse_transform(transform_str::String)

Parse an SVG transform string into a transformation matrix.
"""
function parse_transform(transform_str::String)
    # Start with identity matrix
    matrix = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    
    # Extract individual transform functions
    transform_functions = extract_transform_functions(transform_str)
    
    # Process each transform function in reverse order
    # SVG transformations are applied from right to left
    for transform_func in reverse(transform_functions)
        # Apply the transformation to the matrix
        matrix = apply_transform_function(transform_func, matrix)
    end
    
    return matrix
end

"""
    apply_transform_function(transform_func::String, matrix::Matrix{Float64})

Apply a single transform function to a transformation matrix.
"""
function apply_transform_function(transform_func::String, matrix::Matrix{Float64})
    # Regular expressions to match transform functions
    translate_regex = r"translate\(\s*([^,\)]+)(?:,\s*([^,\)]+))?\s*\)"
    scale_regex = r"scale\(\s*([^,\)]+)(?:,\s*([^,\)]+))?\s*\)"
    rotate_regex = r"rotate\(\s*([^,\)]+)(?:,\s*([^,\)]+)(?:,\s*([^,\)]+))?)?\s*\)"
    matrix_regex = r"matrix\(\s*([^,\s]+)[\s,]+([^,\s]+)[\s,]+([^,\s]+)[\s,]+([^,\s]+)[\s,]+([^,\s]+)[\s,]+([^,\s]+)\s*\)"
    skewX_regex = r"skewX\(\s*([^,\)]+)\s*\)"
    skewY_regex = r"skewY\(\s*([^,\)]+)\s*\)"
    
    # Parse translate
    m = match(translate_regex, transform_func)
    if !isnothing(m)
        tx = parse(Float64, m.captures[1])
        ty = isnothing(m.captures[2]) ? 0.0 : parse(Float64, m.captures[2])
        translate_matrix = [1.0 0.0 tx; 0.0 1.0 ty; 0.0 0.0 1.0]
        return translate_matrix * matrix
    end
    
    # Parse scale
    m = match(scale_regex, transform_func)
    if !isnothing(m)
        sx = parse(Float64, m.captures[1])
        sy = isnothing(m.captures[2]) ? sx : parse(Float64, m.captures[2])
        scale_matrix = [sx 0.0 0.0; 0.0 sy 0.0; 0.0 0.0 1.0]
        return scale_matrix * matrix
    end
    
    # Parse rotate
    m = match(rotate_regex, transform_func)
    if !isnothing(m)
        angle_deg = parse(Float64, m.captures[1])
        angle_rad = deg2rad(angle_deg)
        cos_a = cos(angle_rad)
        sin_a = sin(angle_rad)
        
        if isnothing(m.captures[2]) || isnothing(m.captures[3])
            # Rotate around origin
            rotate_matrix = [cos_a -sin_a 0.0; sin_a cos_a 0.0; 0.0 0.0 1.0]
            return rotate_matrix * matrix
        else
            # Rotate around point (cx, cy)
            cx = parse(Float64, m.captures[2])
            cy = parse(Float64, m.captures[3])
            
            # Translate to origin, rotate, translate back
            translate_to_origin = [1.0 0.0 -cx; 0.0 1.0 -cy; 0.0 0.0 1.0]
            rotate_matrix = [cos_a -sin_a 0.0; sin_a cos_a 0.0; 0.0 0.0 1.0]
            translate_back = [1.0 0.0 cx; 0.0 1.0 cy; 0.0 0.0 1.0]
            
            transform_matrix = translate_back * rotate_matrix * translate_to_origin
            return transform_matrix * matrix
        end
    end
    
    # Parse matrix
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
    
    # Parse skewX
    m = match(skewX_regex, transform_func)
    if !isnothing(m)
        angle_deg = parse(Float64, m.captures[1])
        angle_rad = deg2rad(angle_deg)
        tan_a = tan(angle_rad)
        
        skew_matrix = [1.0 tan_a 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        return skew_matrix * matrix
    end
    
    # Parse skewY
    m = match(skewY_regex, transform_func)
    if !isnothing(m)
        angle_deg = parse(Float64, m.captures[1])
        angle_rad = deg2rad(angle_deg)
        tan_a = tan(angle_rad)
        
        skew_matrix = [1.0 0.0 0.0; tan_a 1.0 0.0; 0.0 0.0 1.0]
        return skew_matrix * matrix
    end
    
    # Unknown transform function
    return matrix
end

"""
    extract_transform_functions(transform_str::String)

Extract individual transform functions from a transform string.
"""
function extract_transform_functions(transform_str::String)
    transform_functions = String[]
    
    # Regular expressions to match transform functions
    transform_regex = r"([a-zA-Z]+\([^)]*\))"
    
    # Find all transform functions in the string
    for m in eachmatch(transform_regex, transform_str)
        push!(transform_functions, m.captures[1])
    end
    
    return transform_functions
end

"""
    matrix_to_transform_str(matrix::Matrix{Float64})

Convert a transformation matrix to an SVG transform string.
"""
function matrix_to_transform_str(matrix::Matrix{Float64})
    # Extract the 6 values from the matrix (a, b, c, d, e, f)
    a = matrix[1, 1]
    b = matrix[2, 1]
    c = matrix[1, 2]
    d = matrix[2, 2]
    e = matrix[1, 3]
    f = matrix[2, 3]
    
    # Format with 6 decimal places to avoid floating point issues
    return @sprintf("matrix(%.6f,%.6f,%.6f,%.6f,%.6f,%.6f)", a, b, c, d, e, f)
end

"""
    shape_to_polyline(shape::XML.Node, num_points::Int)

Convert an SVG shape to a polyline.
"""
function shape_to_polyline(shape::XML.Node, num_points::Int)
    tag = XML.tag(shape)
    attrs = get_attributes(shape)
    
    # Extract style attributes
    style = extract_style_attributes(attrs)
    
    # Convert based on shape type
    if tag == "rect"
        return rect_to_polyline(attrs, style)
    elseif tag == "circle"
        return circle_to_polyline(attrs, style, num_points)
    elseif tag == "ellipse"
        return ellipse_to_polyline(attrs, style, num_points)
    elseif tag == "line"
        return line_to_polyline(attrs, style)
    elseif tag == "polyline"
        return parse_polyline(attrs, style, false)
    elseif tag == "polygon"
        return parse_polyline(attrs, style, true)
    elseif tag == "path"
        return path_to_polyline(attrs, style, num_points)
    else
        error("Unsupported shape type: $tag")
    end
end

"""
    extract_style_attributes(attrs::Dict{String, String})

Extract style attributes from element attributes and convert colors to RGBA arrays.
"""
function extract_style_attributes(attrs::Dict{String, String})
    # Extract basic style properties
    fill_color = get(attrs, "fill", "black")
    stroke_color = get(attrs, "stroke", "none")
    stroke_width_str = get(attrs, "stroke-width", "1")
    opacity_str = get(attrs, "opacity", "1")
    fill_opacity_str = get(attrs, "fill-opacity", "1")
    stroke_opacity_str = get(attrs, "stroke-opacity", "1")
    
    # Parse colors to RGBA arrays
    fill_rgba = parse_color(fill_color)
    stroke_rgba = parse_color(stroke_color)
    
    # Parse opacity values
    opacity = parse(Float64, opacity_str)
    fill_opacity = parse(Float64, fill_opacity_str)
    stroke_opacity = parse(Float64, stroke_opacity_str)
    
    # Apply opacity to colors
    fill_rgba[4] *= fill_opacity * opacity
    stroke_rgba[4] *= stroke_opacity * opacity
    
    # Parse stroke width
    stroke_width = parse_stroke_width(stroke_width_str)
    
    # Collect other style attributes
    other_style = Dict{String, String}()
    for (key, value) in attrs
        if key != "fill" && key != "stroke" && key != "stroke-width" && 
           key != "opacity" && key != "fill-opacity" && key != "stroke-opacity"
            other_style[key] = value
        end
    end
    
    # Create Style struct
    return Style(fill_rgba, stroke_rgba, stroke_width, opacity, other_style)
end

"""
    rect_to_polyline(attrs::Dict{String, String}, style::Style)

Convert a rectangle to a polyline.
"""
function rect_to_polyline(attrs::Dict{String, String}, style::Style)
    # Extract rectangle attributes
    x = parse(Float64, get(attrs, "x", "0"))
    y = parse(Float64, get(attrs, "y", "0"))
    width = parse(Float64, attrs["width"])
    height = parse(Float64, attrs["height"])
    
    # Create points for the rectangle (4 corners)
    points = [
        SVGPoint(x, y),
        SVGPoint(x + width, y),
        SVGPoint(x + width, y + height),
        SVGPoint(x, y + height)
    ]
    
    # Return as a closed polyline (polygon)
    return SVGPoints(points, true, style)
end

"""
    circle_to_polyline(attrs::Dict{String, String}, style::Style, num_points::Int)

Convert a circle to a polyline.
"""
function circle_to_polyline(attrs::Dict{String, String}, style::Style, num_points::Int)
    # Extract circle attributes
    cx = parse(Float64, get(attrs, "cx", "0"))
    cy = parse(Float64, get(attrs, "cy", "0"))
    r = parse(Float64, attrs["r"])
    
    # Create points around the circle
    points = Vector{SVGPoint}(undef, num_points)
    for i in 1:num_points
        angle = 2π * (i - 1) / num_points
        x = cx + r * cos(angle)
        y = cy + r * sin(angle)
        points[i] = SVGPoint(x, y)
    end
    
    # Return as a closed polyline (polygon)
    return SVGPoints(points, true, style)
end

"""
    ellipse_to_polyline(attrs::Dict{String, String}, style::Style, num_points::Int)

Convert an ellipse to a polyline.
"""
function ellipse_to_polyline(attrs::Dict{String, String}, style::Style, num_points::Int)
    # Extract ellipse attributes
    cx = parse(Float64, get(attrs, "cx", "0"))
    cy = parse(Float64, get(attrs, "cy", "0"))
    rx = parse(Float64, attrs["rx"])
    ry = parse(Float64, attrs["ry"])
    
    # Create points around the ellipse
    points = Vector{SVGPoint}(undef, num_points)
    for i in 1:num_points
        angle = 2π * (i - 1) / num_points
        x = cx + rx * cos(angle)
        y = cy + ry * sin(angle)
        points[i] = SVGPoint(x, y)
    end
    
    # Return as a closed polyline (polygon)
    return SVGPoints(points, true, style)
end

"""
    line_to_polyline(attrs::Dict{String, String}, style::Style)

Convert a line to a polyline.
"""
function line_to_polyline(attrs::Dict{String, String}, style::Style)
    # Extract line attributes
    x1 = parse(Float64, get(attrs, "x1", "0"))
    y1 = parse(Float64, get(attrs, "y1", "0"))
    x2 = parse(Float64, get(attrs, "x2", "0"))
    y2 = parse(Float64, get(attrs, "y2", "0"))
    
    # Create points for the line (just 2 points)
    points = [SVGPoint(x1, y1), SVGPoint(x2, y2)]
    
    # Return as an open polyline
    return SVGPoints(points, false, style)
end

"""
    parse_polyline(attrs::Dict{String, String}, style::Style, closed::Bool)

Parse a polyline or polygon from attributes.
"""
function parse_polyline(attrs::Dict{String, String}, style::Style, closed::Bool)
    # Extract points string
    points_str = attrs["points"]
    
    # Parse points
    points = SVGPoint[]
    for point_str in Base.split(points_str)
        if isempty(point_str)
            continue
        end
        
        coords = Base.split(point_str, ',')
        if length(coords) >= 2
            x = parse(Float64, coords[1])
            y = parse(Float64, coords[2])
            push!(points, SVGPoint(x, y))
        end
    end
    
    # Return as a polyline
    return SVGPoints(points, closed, style)
end

"""
    path_to_polyline(attrs::Dict{String, String}, style::Style, num_points::Int)

Convert a path to a polyline.
"""
function path_to_polyline(attrs::Dict{String, String}, style::Style, num_points::Int)
    # Extract path data
    path_data = attrs["d"]
    
    # Parse path data and convert to polylines
    return parse_path_data(path_data, style, num_points)
end

"""
    parse_path_data(path_data::String, style::Style, num_points::Int)

Parse SVG path data and convert to a polyline.
"""
function parse_path_data(path_data::String, style::Style, num_points::Int)
    # Tokenize path data
    tokens = tokenize_path_data(path_data)
    
    # Parse tokens
    points = SVGPoint[]
    current_point = SVGPoint(0.0, 0.0)
    start_point = SVGPoint(0.0, 0.0)
    closed = false
    
    i = 1
    while i <= length(tokens)
        command = tokens[i]
        i += 1
        
        if command == "M" || command == "m"
            # Move to
            relative = command == "m"
            x = parse(Float64, tokens[i])
            y = parse(Float64, tokens[i+1])
            i += 2
            
            if relative
                current_point = SVGPoint(current_point.x + x, current_point.y + y)
            else
                current_point = SVGPoint(x, y)
            end
            
            start_point = current_point
            push!(points, current_point)
            
        elseif command == "L" || command == "l"
            # Line to
            relative = command == "l"
            x = parse(Float64, tokens[i])
            y = parse(Float64, tokens[i+1])
            i += 2
            
            if relative
                current_point = SVGPoint(current_point.x + x, current_point.y + y)
            else
                current_point = SVGPoint(x, y)
            end
            
            push!(points, current_point)
            
        elseif command == "H" || command == "h"
            # Horizontal line to
            relative = command == "h"
            x = parse(Float64, tokens[i])
            i += 1
            
            if relative
                current_point = SVGPoint(current_point.x + x, current_point.y)
            else
                current_point = SVGPoint(x, current_point.y)
            end
            
            push!(points, current_point)
            
        elseif command == "V" || command == "v"
            # Vertical line to
            relative = command == "v"
            y = parse(Float64, tokens[i])
            i += 1
            
            if relative
                current_point = SVGPoint(current_point.x, current_point.y + y)
            else
                current_point = SVGPoint(current_point.x, y)
            end
            
            push!(points, current_point)
            
        elseif command == "C" || command == "c"
            # Cubic Bezier curve
            relative = command == "c"
            
            x1 = parse(Float64, tokens[i])
            y1 = parse(Float64, tokens[i+1])
            x2 = parse(Float64, tokens[i+2])
            y2 = parse(Float64, tokens[i+3])
            x = parse(Float64, tokens[i+4])
            y = parse(Float64, tokens[i+5])
            i += 6
            
            if relative
                x1 += current_point.x
                y1 += current_point.y
                x2 += current_point.x
                y2 += current_point.y
                x += current_point.x
                y += current_point.y
            end
            
            # Sample points along the cubic Bezier curve
            for t in range(0, 1, length=num_points)
                # Cubic Bezier formula
                bx = (1-t)^3 * current_point.x + 3*(1-t)^2*t * x1 + 3*(1-t)*t^2 * x2 + t^3 * x
                by = (1-t)^3 * current_point.y + 3*(1-t)^2*t * y1 + 3*(1-t)*t^2 * y2 + t^3 * y
                
                push!(points, SVGPoint(bx, by))
            end
            
            current_point = SVGPoint(x, y)
            
        elseif command == "S" || command == "s"
            # Smooth cubic Bezier curve
            relative = command == "s"
            
            # Reflect previous control point
            if length(points) >= 2
                prev_control = points[end-1]
                x1 = 2 * current_point.x - prev_control.x
                y1 = 2 * current_point.y - prev_control.y
            else
                x1 = current_point.x
                y1 = current_point.y
            end
            
            x2 = parse(Float64, tokens[i])
            y2 = parse(Float64, tokens[i+1])
            x = parse(Float64, tokens[i+2])
            y = parse(Float64, tokens[i+3])
            i += 4
            
            if relative
                x2 += current_point.x
                y2 += current_point.y
                x += current_point.x
                y += current_point.y
            end
            
            # Sample points along the cubic Bezier curve
            for t in range(0, 1, length=num_points)
                # Cubic Bezier formula
                bx = (1-t)^3 * current_point.x + 3*(1-t)^2*t * x1 + 3*(1-t)*t^2 * x2 + t^3 * x
                by = (1-t)^3 * current_point.y + 3*(1-t)^2*t * y1 + 3*(1-t)*t^2 * y2 + t^3 * y
                
                push!(points, SVGPoint(bx, by))
            end
            
            current_point = SVGPoint(x, y)
            
        elseif command == "Q" || command == "q"
            # Quadratic Bezier curve
            relative = command == "q"
            
            x1 = parse(Float64, tokens[i])
            y1 = parse(Float64, tokens[i+1])
            x = parse(Float64, tokens[i+2])
            y = parse(Float64, tokens[i+3])
            i += 4
            
            if relative
                x1 += current_point.x
                y1 += current_point.y
                x += current_point.x
                y += current_point.y
            end
            
            # Sample points along the quadratic Bezier curve
            for t in range(0, 1, length=num_points)
                # Quadratic Bezier formula
                bx = (1-t)^2 * current_point.x + 2*(1-t)*t * x1 + t^2 * x
                by = (1-t)^2 * current_point.y + 2*(1-t)*t * y1 + t^2 * y
                
                push!(points, SVGPoint(bx, by))
            end
            
            current_point = SVGPoint(x, y)
            
        elseif command == "T" || command == "t"
            # Smooth quadratic Bezier curve
            relative = command == "t"
            
            # Reflect previous control point
            if length(points) >= 2
                prev_control = points[end-1]
                x1 = 2 * current_point.x - prev_control.x
                y1 = 2 * current_point.y - prev_control.y
            else
                x1 = current_point.x
                y1 = current_point.y
            end
            
            x = parse(Float64, tokens[i])
            y = parse(Float64, tokens[i+1])
            i += 2
            
            if relative
                x += current_point.x
                y += current_point.y
            end
            
            # Sample points along the quadratic Bezier curve
            for t in range(0, 1, length=num_points)
                # Quadratic Bezier formula
                bx = (1-t)^2 * current_point.x + 2*(1-t)*t * x1 + t^2 * x
                by = (1-t)^2 * current_point.y + 2*(1-t)*t * y1 + t^2 * y
                
                push!(points, SVGPoint(bx, by))
            end
            
            current_point = SVGPoint(x, y)
            
        elseif command == "A" || command == "a"
            # Elliptical arc
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
                x += current_point.x
                y += current_point.y
            end
            
            # Convert elliptical arc to polyline
            arc_points = elliptical_arc_to_points(
                current_point.x, current_point.y,
                rx, ry,
                deg2rad(x_axis_rotation),
                large_arc_flag != 0,
                sweep_flag != 0,
                x, y,
                num_points
            )
            
            append!(points, arc_points)
            current_point = SVGPoint(x, y)
            
        elseif command == "Z" || command == "z"
            # Close path
            if length(points) > 0
                push!(points, start_point)
                closed = true
            end
            
        else
            # Unknown command
            @warn "Unknown path command: $command"
        end
    end
    
    return SVGPoints(points, closed, style)
end

"""
    tokenize_path_data(path_data::String)

Tokenize SVG path data into commands and parameters.
"""
function tokenize_path_data(path_data::String)
    # Replace commas with spaces
    path_data = replace(path_data, "," => " ")
    
    # Add spaces around commands
    for cmd in "MmLlHhVvCcSsQqTtAaZz"
        path_data = replace(path_data, cmd => " $cmd ")
    end
    
    # Split into tokens and filter out empty strings
    tokens = filter(!isempty, Base.split(path_data))
    
    return tokens
end

"""
    elliptical_arc_to_points(x1, y1, rx, ry, phi, large_arc, sweep, x2, y2, num_points)

Convert an elliptical arc to a sequence of points.
"""
function elliptical_arc_to_points(x1, y1, rx, ry, phi, large_arc, sweep, x2, y2, num_points)
    # Implementation based on SVG spec: https://www.w3.org/TR/SVG/implnote.html#ArcImplementationNotes
    
    # Ensure radii are positive
    rx = abs(rx)
    ry = abs(ry)
    
    # If radii are zero, treat as a straight line
    if rx < 1e-10 || ry < 1e-10
        return [SVGPoint(x2, y2)]
    end
    
    # Step 1: Transform to origin
    cos_phi = cos(phi)
    sin_phi = sin(phi)
    
    # Calculate midpoint
    dx = (x1 - x2) / 2
    dy = (y1 - y2) / 2
    
    # Transform to origin
    x1p = cos_phi * dx + sin_phi * dy
    y1p = -sin_phi * dx + cos_phi * dy
    
    # Ensure radii are large enough
    lambda = (x1p / rx)^2 + (y1p / ry)^2
    if lambda > 1
        rx *= sqrt(lambda)
        ry *= sqrt(lambda)
    end
    
    # Step 2: Calculate center
    sign = large_arc == sweep ? -1 : 1
    sq = ((rx*ry)^2 - (rx*y1p)^2 - (ry*x1p)^2) / ((rx*y1p)^2 + (ry*x1p)^2)
    sq = max(0, sq)  # Ensure non-negative
    coef = sign * sqrt(sq)
    
    cxp = coef * (rx * y1p) / ry
    cyp = coef * (-ry * x1p) / rx
    
    # Step 3: Transform back
    cx = cos_phi * cxp - sin_phi * cyp + (x1 + x2) / 2
    cy = sin_phi * cxp + cos_phi * cyp + (y1 + y2) / 2
    
    # Step 4: Calculate angles
    ux = (x1p - cxp) / rx
    uy = (y1p - cyp) / ry
    vx = (-x1p - cxp) / rx
    vy = (-y1p - cyp) / ry
    
    # Calculate start angle
    theta1 = angle_between(1.0, 0.0, ux, uy)
    
    # Calculate angle delta
    delta = angle_between(ux, uy, vx, vy)
    if !sweep && delta > 0
        delta -= 2π
    elseif sweep && delta < 0
        delta += 2π
    end
    
    # Generate points along the arc
    points = SVGPoint[]
    for i in 0:num_points
        t = theta1 + delta * i / num_points
        
        # Calculate point on unit circle
        px = cos(t)
        py = sin(t)
        
        # Transform to ellipse
        x = cos_phi * rx * px - sin_phi * ry * py + cx
        y = sin_phi * rx * px + cos_phi * ry * py + cy
        
        push!(points, SVGPoint(x, y))
    end
    
    return points
end

"""
    angle_between(ux, uy, vx, vy)

Calculate the angle between two vectors.
"""
function angle_between(ux, uy, vx, vy)
    # Calculate dot product
    dot = ux * vx + uy * vy
    
    # Calculate magnitudes
    len_u = sqrt(ux^2 + uy^2)
    len_v = sqrt(vx^2 + vy^2)
    
    # Calculate angle
    angle = acos(clamp(dot / (len_u * len_v), -1.0, 1.0))
    
    # Determine sign
    if ux * vy - uy * vx < 0
        angle = -angle
    end
    
    return angle
end

"""
    style_to_dict(style::Style)

Convert a Style struct to a dictionary for JSON serialization.
"""
function style_to_dict(style::Style)
    # Create a dictionary with all style properties
    style_dict = Dict{String, Any}(
        "fill" => style.fill,
        "stroke" => style.stroke,
        "stroke_width" => style.stroke_width,
        "opacity" => style.opacity
    )
    
    # Add other style properties
    for (key, value) in style.other
        style_dict[key] = value
    end
    
    return style_dict
end

"""
    node_to_dict(node::SVGNode)

Convert an SVG node to a dictionary for JSON serialization.
"""
function node_to_dict(node::SVGNode)
    if node isa SVGGroup
        # Convert group node
        children = [node_to_dict(child) for child in node.children]
        return Dict(
            "type" => "group",
            "id" => node.id,
            "children" => children
        )
    elseif node isa SVGPolyline
        # Convert polyline node
        polyline = node.polyline
        points = [[p.x, p.y] for p in polyline.points]
        return Dict(
            "type" => "polyline",
            "id" => node.id,
            "closed" => polyline.closed,
            "style" => style_to_dict(polyline.style),
            "points" => points
        )
    else
        # Unknown node type
        return Dict("type" => "unknown")
    end
end

"""
    save_as_json(node::SVGNode, filename::String)

Save an SVG node tree as a JSON file.
"""
function save_as_json(node::SVGNode, filename::String)
    # Convert to dictionary
    dict = node_to_dict(node)
    
    # Write to file
    open(filename, "w") do f
        write(f, dict_to_json(dict))
    end
    
    return filename
end

"""
    dict_to_json(d, indent=0)

Convert a dictionary to a JSON string.
"""
function dict_to_json(d, indent=0)
    indent_str = " " ^ indent
    next_indent_str = " " ^ (indent + 2)
    
    if d isa Dict
        # Convert dictionary
        parts = String[]
        push!(parts, "{\n")
        
        for (i, (key, value)) in enumerate(d)
            push!(parts, next_indent_str * "\"$(key)\": ")
            push!(parts, dict_to_json(value, indent + 2))
            
            if i < length(d)
                push!(parts, ",\n")
            else
                push!(parts, "\n")
            end
        end
        
        push!(parts, indent_str * "}")
        return join(parts, "")
    elseif d isa Vector
        # Convert array
        if isempty(d)
            return "[]"
        end
        
        parts = String[]
        push!(parts, "[\n")
        
        for (i, value) in enumerate(d)
            push!(parts, next_indent_str)
            push!(parts, dict_to_json(value, indent + 2))
            
            if i < length(d)
                push!(parts, ",\n")
            else
                push!(parts, "\n")
            end
        end
        
        push!(parts, indent_str * "]")
        return join(parts, "")
    elseif d isa String
        # Convert string
        return "\"$(escape_json_string(d))\""
    elseif d isa Bool || d isa Number || d === nothing
        # Convert primitive
        return string(d)
    else
        # Convert other types
        return "\"$(string(d))\""
    end
end

"""
    escape_json_string(s)

Escape special characters in a string for JSON.
"""
function escape_json_string(s)
    replacements = [
        ('\\', "\\\\"),
        ('"', "\\\""),
        ('\b', "\\b"),
        ('\f', "\\f"),
        ('\n', "\\n"),
        ('\r', "\\r"),
        ('\t', "\\t")
    ]
    
    result = s
    for (char, replacement) in replacements
        result = replace(result, char => replacement)
    end
    
    return result
end



# ////////////////////////////////////////////////////////////////
# for ALBERTO
export svg_to_plasm
function svg_to_plasm(node)

    if node isa SVGGroup
        v=[]
        for child in node.children
          push!(v, svg_to_plasm(child))
        end
        return length(v) > 0 ? STRUCT(v) : Hpc()
    
    elseif node isa SVGPolyline
        polyline = node.polyline
        
        line_color=Point4d(polyline.style.stroke)
        line_width=max(1,Int(round(polyline.style.stroke_width)))
        face_color=Point4d(polyline.style.fill)
  
        # println(line_color, line_width, face_color)
  
        points=[[p.x,p.y] for p in polyline.points]
        if polyline.closed
  
          # this is in general wrong (could be non-convex) 
          # ret=MKPOL(points, [collect(1:length(points))], [[1]]) 
          
          push!(points, points[1]) # close the polygon
          ret=MKPOL(points, [[i, i + 1] for i in 1:length(points)-1])
        else
          ret=MKPOL(points, [[i, i + 1] for i in 1:length(points)-1])
        end
        
        return PROPERTIES(ret, 
          Properties(
            "line_color" => line_color,
            "line_width" => line_width,
            "face_color" => face_color
          ))
    end
  end
  
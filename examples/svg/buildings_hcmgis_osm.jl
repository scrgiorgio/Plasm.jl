
using Plasm
import Base:*
using CSV
using DataFrames

svg_node=parse_svg_file("./examples/svg/buildings_hcmgis_osm.svg", num_points=24)
obj=svg_to_plasm(svg_node)
VIEW(obj)
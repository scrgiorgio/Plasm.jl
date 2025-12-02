using Plasm

#=
svg_to_plasm(node; triangulate=nothing, output_type=:hpc)

triangulate:
  - nothing: edges only (wireframe)
  - :global: triangulate all shapes together (handles holes/overlaps)
  - :local: triangulate each shape independently (faster, simpler shapes only)

output_type:
  - :hpc: Hierarchical Polyhedral Complex
  - :lar: Linear Algebraic Representation
=#

holes = parse_svg_file("./test/svg/holes.svg")
buildings = parse_svg_file("./test/svg/buildings.svg")

# Hpc output
begin

  # just wireframe
  VIEW(svg_to_plasm(parse_svg_file("./test/svg/apartment.svg")))
  VIEW(svg_to_plasm(parse_svg_file("./test/svg/curved-1.svg")))
  VIEW(svg_to_plasm(parse_svg_file("./test/svg/curved-2.svg")))
  VIEW(svg_to_plasm(parse_svg_file("./test/svg/curved-3.svg")))
  VIEW(svg_to_plasm(parse_svg_file("./test/svg/curved-4.svg")))
  VIEW(svg_to_plasm(parse_svg_file("./test/svg/example.svg")))
  VIEW(svg_to_plasm(parse_svg_file("./test/svg/test-1.svg")))
  VIEW(svg_to_plasm(parse_svg_file("./test/svg/test-2.svg")))
  VIEW(svg_to_plasm(parse_svg_file("./test/svg/test-3.svg")))
  VIEW(svg_to_plasm(parse_svg_file("./test/svg/test-4.svg")))
  VIEW(svg_to_plasm(parse_svg_file("./test/svg/test-5.svg")))
  VIEW(svg_to_plasm(parse_svg_file("./test/svg/test-6.svg")))
  VIEW(svg_to_plasm(parse_svg_file("./test/svg/test-7.svg")))

  
  # === global triangulation ===
  hpc::Hpc = svg_to_plasm(holes, triangulate=:global, output_type=:hpc)
  VIEW(hpc, title="HPC from global triangulation")

  # === local triangulation (this is possible only if you are sure each building is well defined) ===
  hpc::Hpc = svg_to_plasm(buildings, triangulate=:local, output_type=:hpc)
  VIEW(hpc, title="HPC from local triangulation")
end

# LAR output
begin

  # just wireframe
  VIEWCOMPLEX(svg_to_plasm(parse_svg_file("./test/svg/apartment.svg") , output_type=:lar))

  # === global triangulation ===
  lar::Lar = svg_to_plasm(holes, triangulate=:global, output_type=:lar)
  VIEWCOMPLEX(lar, title="LAR from global triangulation")

  # === local triangulation (this is possible only if you are sure each building is well defined) ===
  lar::Lar = svg_to_plasm(buildings, triangulate=:local, output_type=:lar)
  VIEWCOMPLEX(lar, title="LAR from local triangulation")
end












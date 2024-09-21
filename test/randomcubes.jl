using Plasm
using Random

Random.seed!(0)

hpc = STRUCT([RandomCube(0.2,2.0) for I in 1:6])
lar = LAR(hpc)

arrangement = ARRANGE3D(lar)

arrangement = SIMPLIFY(arrangement)

arrangement=without_outer_atom(arrangement)

# if you want to see the atoms...
show_atoms=false
if show_atoms
  for (A, sel) in enumerate(arrangement.C[:CF])
    atom = SELECT(arrangement, sel)
    @show(atom)
    VIEWCOMPLEX(atom, explode=[1.4, 1.4, 1.4], show=["V", "EV", "FV", "V_text", "EV_text", "FV_text"], face_color=TRANSPARENT)
  end
end

# show faces, exploding each face by its centroid
VIEWCOMPLEX(arrangement, show=["FV"], explode=[1.4,1.4,1.4])

# show faces, but keep the atom together
VIEWCOMPLEX(arrangement, show=["FV","atom"], explode=[1.4,1.4,1.4])
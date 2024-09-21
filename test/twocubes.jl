using Plasm
using Random

Random.seed!(0)

cube = CUBOID([1,1,1])
hpc = STRUCT(
  cube, 
  T(1,2,3)(.5,.5,0.75), 
  R(2,3)(pi/4), 
  R(1,3)(pi/4),
  cube)
lar = LAR(hpc)

arrangement = ARRANGE3D(lar)

arrangement = SIMPLIFY(arrangement)

# do not want to see outer atom
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

# show faces, exploding each face by sits centroid
VIEWCOMPLEX(arrangement, show=["FV"], explode=[1.2,1.2,2.0])

# show faces, but keep the atom together
VIEWCOMPLEX(arrangement, show=["FV","atom"], explode=[1.2,1.2,2.0])
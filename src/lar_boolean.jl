
export settestpoints, testinternalpoint, getinternalpoint, internalpoints, rayintersection, planemap, bool3d

################################################################################
#	After the arrangement, extract all the d-cells from (d-1)-coboundary as isolated polyhedra.
#	Then compute a single interior point for each of them.
#	Then compare each such point against all input boundaries, in order to compute those which it was interior to. Extend this point membership as 3-cell containment within the relative input solids.
#	The point membership with a boundary consists in the parity count of the intersection points of a vertical ray starting at the test point, with the boundary surface.
#
################################################################################

# //////////////////////////////////////////////////////////////////////////////
# working in 2D ///////////////////////////////////////////////////////////////
export pointInPolygonClassification
""" Test di contenimento di un punto 2D in un poligono 2D """
function pointInPolygonClassification(V, EV)

  function setTile(box)
    tiles = [[9, 1, 5], [8, 0, 4], [10, 2, 6]]
    b1, b2, b3, b4 = box
    """ code point position w.r.t query box using Bitwise OR """
    function tileCode(point)
      x, y = point
      code = 0
      if y > b1
        code = code | 1
      end
      if y < b2
        code = code | 2
      end
      if x > b3
        code = code | 4
      end
      if x < b4
        code = code | 8
      end
      return code
    end
    return tileCode
  end

  function pointInPolygonClassification0(pnt)
    x, y = pnt
    xmin, xmax, ymin, ymax = x, x, y, y
    tilecode = setTile([ymax, ymin, xmax, xmin])
    count, status = 0, 0

    for (k, edge) in enumerate(EV)
      p1, p2 = V[:, edge[1]], V[:, edge[2]]
      (x1, y1), (x2, y2) = p1, p2
      c1, c2 = tilecode(p1), tilecode(p2)
      c_edge, c_un, c_int = xor(c1, c2), c1 | c2, c1 & c2

      if (c_edge == 0) & (c_un == 0)
        return "p_on"
      elseif (c_edge == 12) & (c_un == c_edge)
        return "p_on"
      elseif c_edge == 3
        if c_int == 0
          return "p_on"
        elseif c_int == 4
          count += 1
        end
      elseif c_edge == 15
        x_int = ((y - y2) * (x1 - x2) / (y1 - y2)) + x2
        if x_int > x
          count += 1
        elseif x_int == x
          return "p_on"
        end
      elseif (c_edge == 13) & ((c1 == 4) | (c2 == 4))
        crossingTest(1, 2, status, count)
      elseif (c_edge == 14) & ((c1 == 4) | (c2 == 4))
        crossingTest(2, 1, status, count)
      elseif c_edge == 7
        count += 1
      elseif c_edge == 11
        count = count
      elseif c_edge == 1
        if c_int == 0
          return "p_on"
        elseif c_int == 4
          crossingTest(1, 2, status, count)
        end
      elseif c_edge == 2
        if c_int == 0
          return "p_on"
        elseif c_int == 4
          crossingTest(2, 1, status, count)
        end
      elseif (c_edge == 4) & (c_un == c_edge)
        return "p_on"
      elseif (c_edge == 8) & (c_un == c_edge)
        return "p_on"
      elseif c_edge == 5
        if (c1 == 0) | (c2 == 0)
          return "p_on"
        else
          crossingTest(1, 2, status, count)
        end
      elseif c_edge == 6
        if (c1 == 0) | (c2 == 0)
          return "p_on"
        else
          crossingTest(2, 1, status, count)
        end
      elseif (c_edge == 9) & ((c1 == 0) | (c2 == 0))
        return "p_on"
      elseif (c_edge == 10) & ((c1 == 0) | (c2 == 0))
        return "p_on"
      end
    end
    if (round(count) % 2) == 1
      return "p_in"
    else
      return "p_out"
    end
  end
  return pointInPolygonClassification0
end

# //// look for a pair of test point for given atom ////////////////////////////
# working in 3D ////////////////////////////////////////////////////////////////
function settestpoints(V, EV, FE, FV, Fs, copEV, copFE) # V by row
  f = Fs[1]
  e = findnz(copFE[f, :])[1][1] # first (global) edge of first (global) face
  # f,e relative to atom
  f1, f2 = findnz(copFE[f, :])[1][1:2] # two (global) faces incident on it (atom)
  v1, v2 = findnz(copEV[e, :])[1][1:2] # two (global) verts incident on it (atom)
  fdict = Dict(zip(Fs, 1:length(Fs))) # enumerate atom faces (local codes)
  V1 = FV[fdict[f1]]
  V2 = FV[fdict[f2]]
  v1, v2 = intersect(V1, V2) # verified ... !
  # two tringles sharing a common edge 
  t1 = V[:, v1], V[:, v2], V[:, [v for v in V1 if v ≠ v1 && v ≠ v2][1]]# t1 ⊆ ∂atom ??
  t2 = V[:, v2], V[:, v1], V[:, [v for v in V2 if v ≠ v1 && v ≠ v2][1]]# t2 ⊆ ∂atom ??
  n1 = LinearAlgebra.normalize(cross(t1[2] - t1[1], t1[3] - t1[1]))
  n2 = LinearAlgebra.normalize(cross(t2[2] - t2[1], t2[3] - t2[1]))
  p0 = (V[:, v1] + V[:, v2]) ./ 2 # mean point
  n = n1 + n2  # mean normal
  ϵ = 1.0e-6  # Alberto 2024-03-03
  #ϵ = 1.0e-8
  ptest1 = p0 + ϵ * n
  ptest2 = p0 - ϵ * n
  # it is not known if angle(f1,f2) < π;  hence return two opposite points
  return ptest1, ptest2
end
# //////////////////////////////////////////////////////////////////////////////
# working in 3D ////////////////////////////////////////////////////////////////


# //////////////////////////////////////////////////////////////////////////////
""" V,EV,FV original faces (before the space arrangement) """
function testinternalpoint(V, EV, FV)

    """ Filter preparation to select input faces with possible interction with vertical test ray """
  function spaceindex2(point3d::Array{Float64,1})::Function
    function spaceindex0(model)::Array{Int,1}
      V, CV = copy(model[1]), copy(model[2])
      V = [V point3d]
      dim, idx = size(V)
      push!(CV, [idx, idx, idx])
      cellpoints = [V[:, CV[k]]::Points for k = 1:length(CV)]
      #----------------------------------------------------------
      bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints] #bound boxes
      xboxdict = coordintervals(1, bboxes)
      yboxdict = coordintervals(2, bboxes)
      # xs,ys are IntervalTree type
      xs = IntervalTrees.IntervalMap{Float64,Array}()
      for (key, boxset) in xboxdict
        xs[tuple(key...)] = boxset
      end
      ys = IntervalTrees.IntervalMap{Float64,Array}()
      for (key, boxset) in yboxdict
        ys[tuple(key...)] = boxset
      end
      xcovers = boxcovering(bboxes, 1, xs)
      ycovers = boxcovering(bboxes, 2, ys)
      covers = [intersect(pair...) for pair in zip(xcovers, ycovers)]

      # add new code part

      # remove each cell from its cover
      pointcover = setdiff(covers[end], [idx + 1])
      return pointcover[1:end-1]
    end
    return spaceindex0
  end
  copEV = lar2cop(EV)
  copFV = lar2cop(FV)
  copFE = copFV * copEV' .÷ 2
  copFE = convert(ChainOp, copFE)
  function testinternalpoint0(testpoint)
    intersectedfaces = Int64[]
    # spatial index for possible intersections with ray
    faces = spaceindex2(testpoint)((V, FV))
    depot = []
    # face in faces :  indices of faces of possible intersection with ray
    for face in faces
      value = rayintersection(testpoint)(V, FV, face)
      if typeof(value) == Array{Float64,1}
        push!(depot, (face, value))
      end
    end
    # actual containment test of ray point in faces within depot
    for (face, point3d) in depot
      vs, edges, point2d = planemap(V, copEV, copFE, face)(point3d)
      classify = pointInPolygonClassification(vs, edges)
      inOut = classify(point2d)
      if inOut != "p_out"
        push!(intersectedfaces, face)
      end
    end
    return intersectedfaces # faces intersected by vertical ray
  end
  return testinternalpoint0
end


# //////////////////////////////////////////////////////////////////////////////

""" return two opposite internal/external points in an atom """
function getinternalpoint(V, EV, FE, FV, Fs, copEV, copFE) # V by rows
  # look at edges for v1=FV[1][1]
  ptest1, ptest2 = settestpoints(V, EV, FE, FV, Fs, copEV, copFE)

  intersectedfaces = Int64[]
  # for each test point compute the face planes intersected by vertical ray
  dep1, dep2 = [], [] # to store the pairs (face, 3D-point)
  # face in Fs : global indices of faces of current solid
  for (f, face) in enumerate(Fs)
    ret1 = rayintersection(ptest1)(V, FV, f)
    ret2 = rayintersection(ptest2)(V, FV, f)
    if typeof(ret1) == Array{Float64,1}
      push!(dep1, (face, ret1))
    end
    if typeof(ret2) == Array{Float64,1}
      push!(dep2, (face, ret2))
    end
  end
  # transform each plane in 2D and look whether the point is internal
  # return the test point with odd numeber of ray intersections
  k1, k2 = 0, 0
  for (face, point3d) in dep1 # testing ptest1 for interior of atom
    vs, edges, point2d = planemap(V, copEV, copFE, face)(point3d)
    classify = pointInPolygonClassification(vs, edges)
    inOut = classify(point2d)
    if inOut != "p_out"
      k1 += 1
      push!(intersectedfaces, face)
    end
  end
  if k1 % 2 == 1
    return ptest1, intersectedfaces
  else # testing ptest2 for interior of atom
    for (face, point3d) in dep2
      vs, edges, point2d = planemap(V, copEV, copFE, face)(point3d)
      classify = pointInPolygonClassification(vs, edges)
      inOut = classify(point2d)
      if inOut != "p_out"
        k2 += 1
        push!(intersectedfaces, face)
      end
    end
    if k2 % 2 == 1
      return ptest2, intersectedfaces
    else
      error("tertium non datur")
    end
  end
end

# //////////////////////////////////////////////////////////////////////////////
export get_atoms
""" Build 3D polyhedrons in pols from arrangement output (outer+atoms) """
function get_atoms(copEV, copFE, copCF)
  FEs = Array{Array{Int64,1},1}[]
  EVs = Array{Array{Array{Int64,1},1},1}[]
  FVs = Array{Array{Int64,1},1}[]

  CF = [findnz(copCF[k, :])[1] for k = 1:copCF.m] # faces per cell
  FE = [findnz(copFE[k, :])[1] for k = 1:copFE.m] # edges per face
  EV = [findnz(copEV[k, :])[1] for k = 1:copEV.m] # vertices per edge
  FV = [union(CAT([EV[e] for e in f])) for f in FE] # vertices per face

  for k = 1:copCF.m
    push!(FEs, [collect(Set([e for e in FE[f]])) for f in CF[k]])
    push!(EVs, [[EV[e] for e in FE[f]] for f in CF[k]])
    push!(FVs, [collect(Set(vcat([EV[e] for e in FE[f]]...))) for f in CF[k]])
  end

  pols = collect(zip(EVs, FVs, FEs)) # all atoms w global numbering
  return pols, CF
end

# //////////////////////////////////////////////////////////////////////////////
function internalpoints(V, copEV, copFE, copCF) 

  # transform each 3-cell in a solid (via model)
  pols, CF = get_atoms(copEV, copFE, copCF) 
  outerspace, atoms = SELECTATOMS(V, pols)

  # compute, for each `pol` (3-cell, i.e atom) in `pols`, one `internalpoint`
  innerpoints = []
  intersectedfaces = []
  for k = 1:length(atoms)
    (EV, FV, FE), Fs = atoms[k], CF[k]
    EV = union(CAT(EV))
    point, facenumber = getinternalpoint(V, EV, FE, FV, Fs, copEV, copFE)
    push!(innerpoints, point)
    push!(intersectedfaces, facenumber)
  end
  return innerpoints, intersectedfaces
end

# //////////////////////////////////////////////////////////////////////////////
function rayintersection(point3d)
  function rayintersection0(V, FV, face::Int)
    l0, l = point3d, [0, 0, 1.0]
    ps = V[:, FV[face]]  # face points
    p0 = ps[:, 1]
    v1, v2 = ps[:, 2] - p0, ps[:, 3] - p0
    n = LinearAlgebra.normalize(cross(v1, v2))

    denom = dot(n, l)
    if (abs(denom) > 1e-8) #1e-6
      p0l0 = p0 - l0
      t = dot(p0l0, n) / denom
      if t > 0
        return l0 + t * l   # cosa ritorna nel caso else TODO 
      end
    else
      # "ray and face are parallel"
      return false
    end
  end
  return rayintersection0
end

# //////////////////////////////////////////////////////////////////////////////
""" Take model,face,point on face plane and return 2Dface, edges, and point """

function planemap(V, copEV, copFE, face)
  fv, edges = vcycle(copEV, copFE, face)
  # Fv = Dict(zip(1:length(fv),fv))
  # edges = [[Fv[u],Fv[v]] for (u,v) in edges]
  function planemap0(point)
    vs = V[:, fv]
    # Plasm.view(Plasm.numbering(0.5)((vs,[[[k] for k=1:4],edges])))
    #translation
    point = point .- vs[:, 1]
    vs = vs .- vs[:, 1]
    u, v = edges[1]
    z, w = [[z, w] for (z, w) in edges if z == v][1]
    v1 = vs[:, u] - vs[:, v]
    v2 = vs[:, w] - vs[:, v]
    v3 = cross(v2, v1)
    M = [v1 v2 v3]
    vs = inv(M) * [vs point]
    outvs = vs[1:2, 1:end-1]
    outpoint = vs[1:2, end]
    return outvs, edges, outpoint
  end
  return planemap0
end


# //////////////////////////////////////////////////////////////////////////////
function bool3d(assembly::Hpc, W, copEV, copFE, copCF) # W by rows
  innerpoints, _ = internalpoints(V, copEV, copFE, copCF)

  # associate internal points to (original) faces of 3-cells
  listOfLar = AA(LAR)(TOPOS(assembly)) # correspond to evalStruct(assembly)

  function __new2old(lar)
    V, FV, EV = lar.V, lar.C[:FV], lar.C[:EV]
  end
  inputfacenumbers = AA(LEN ∘ S2 ∘ __new2old)(listOfLar)
  cumulative = cumsum([0; inputfacenumbers]) .+ 1
  fspans = collect(zip(cumulative[1:end-1], cumulative[2:end] .- 1))
  span(h) = [j for j = 1:length(fspans) if fspans[j][1] <= h <= fspans[j][2]]  # ??

  # test input data for containment of reference points
  lar=LAR(assembly)
  V, FV, EV = lar.V, lar.C[:FV], lar.C[:EV]
  containmenttest = testinternalpoint(V, EV, FV)

  # currently copCF contains the outercell in first column ...
  # TODO remove first row and column, in case (look at file src/)
  boolmatrix = BitArray(undef, length(innerpoints) + 1, length(fspans) + 1)
  boolmatrix[1, 1] = 1
  for (k, thepoint) in enumerate(innerpoints) # k runs on rows
    faces = containmenttest(thepoint) # items of a row
    count = zeros(Int, length(fspans))
    for h in faces
      for s in span(h)
        count[s] += 1
      end
    end
    isinternal = count .% 2
    boolmatrix[k+1, 2:end] = isinternal
  end
  return boolmatrix
end




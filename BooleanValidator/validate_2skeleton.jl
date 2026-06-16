module Validate2Skeleton

using SparseArrays
using LinearAlgebra

export edge_boundary_matrix,
       validate_EV,
       validate_face_cycles,
       solve_face_flip_balance,
       validate_2skeleton

"""
    edge_boundary_matrix(EV, nv)

Build the signed vertex-edge incidence matrix B1 ∈ {-1,0,1}^{nv × ne}
from an oriented edge list EV, where each edge is [u,v].
The convention is:
    B1[u,e] = -1
    B1[v,e] = +1
"""
function edge_boundary_matrix(EV::Vector, nv::Int)
    ne = length(EV)
    I = Int[]
    J = Int[]
    X = Int8[]
    for (e, edge) in enumerate(EV)
        length(edge) == 2 || throw(ArgumentError("Edge $e does not have 2 endpoints"))
        u, v = edge
        push!(I, u); push!(J, e); push!(X, Int8(-1))
        push!(I, v); push!(J, e); push!(X, Int8(+1))
    end
    return sparse(I, J, X, nv, ne)
end

"""
    graph_components(adjacency)

Connected components of an undirected graph given by adjacency lists.
Returns a vector of components, each component being a vector of local vertex ids.
"""
function graph_components(adjacency::Vector{Vector{Int}})
    n = length(adjacency)
    seen = falses(n)
    comps = Vector{Vector{Int}}()
    for s in 1:n
        if seen[s]
            continue
        end
        comp = Int[]
        stack = [s]
        seen[s] = true
        while !isempty(stack)
            u = pop!(stack)
            push!(comp, u)
            for v in adjacency[u]
                if !seen[v]
                    seen[v] = true
                    push!(stack, v)
                end
            end
        end
        push!(comps, sort(comp))
    end
    return comps
end

"""
    validate_EV(V, EV; expect_no_isolated=true)

Basic algebraic/combinatorial validation of the 1-skeleton.
Returns a named tuple with diagnostics and, if valid, the signed matrix B1.
"""
function validate_EV(V::AbstractMatrix, EV::Vector; expect_no_isolated::Bool=true)
    nv = size(V, 2)
    ne = length(EV)

    errors = String[]
    warnings = String[]

    bad_length = Int[]
    degenerate = Int[]
    out_of_range = Int[]

    duplicate_map = Dict{Tuple{Int,Int}, Vector{Int}}()

    for (e, edge) in enumerate(EV)
        if length(edge) != 2
            push!(bad_length, e)
            continue
        end
        u, v = edge
        if !(1 <= u <= nv && 1 <= v <= nv)
            push!(out_of_range, e)
            continue
        end
        if u == v
            push!(degenerate, e)
        end
        key = minmax(u, v)
        if haskey(duplicate_map, key)
            push!(duplicate_map[key], e)
        else
            duplicate_map[key] = [e]
        end
    end

    duplicates = Dict(k => v for (k,v) in duplicate_map if length(v) > 1)

    !isempty(bad_length) && push!(errors, "edges with wrong number of endpoints: $(bad_length)")
    !isempty(out_of_range) && push!(errors, "edges with endpoint indices out of range: $(out_of_range)")
    !isempty(degenerate) && push!(errors, "degenerate edges with coincident endpoints: $(degenerate)")
    !isempty(duplicates) && push!(warnings, "duplicate unordered edges found: $(duplicates)")

    B1 = nothing
    if isempty(bad_length) && isempty(out_of_range)
        B1 = edge_boundary_matrix(EV, nv)

        bad_nnz_cols = Int[]
        bad_sign_cols = Int[]
        bad_sum_cols = Int[]
        for e in 1:ne
            _, vals = findnz(B1[:, e])
            if length(vals) != 2
                push!(bad_nnz_cols, e)
            end
            if sort(collect(vals)) != Int8[-1, +1]
                push!(bad_sign_cols, e)
            end
            if sum(vals) != 0
                push!(bad_sum_cols, e)
            end
        end

        !isempty(bad_nnz_cols) && push!(errors, "columns of B1 with nnz != 2: $(bad_nnz_cols)")
        !isempty(bad_sign_cols) && push!(errors, "columns of B1 without signs {-1,+1}: $(bad_sign_cols)")
        !isempty(bad_sum_cols) && push!(errors, "columns of B1 whose signed sum is not zero: $(bad_sum_cols)")
    end

    degree = zeros(Int, nv)
    adjacency = [Int[] for _ in 1:nv]
    for edge in EV
        if length(edge) == 2
            u, v = edge
            if 1 <= u <= nv && 1 <= v <= nv && u != v
                degree[u] += 1
                degree[v] += 1
                push!(adjacency[u], v)
                push!(adjacency[v], u)
            end
        end
    end

    isolated = findall(==(0), degree)
    dangling = findall(==(1), degree)

    if expect_no_isolated && !isempty(isolated)
        push!(warnings, "isolated vertices: $(isolated)")
    end
    if !isempty(dangling)
        push!(warnings, "degree-1 vertices: $(dangling)")
    end

    components = graph_components(adjacency)
    if length(components) > 1
        push!(warnings, "graph has $(length(components)) connected components: $(components)")
    end

    ok = isempty(errors)

    return (
        ok = ok,
        errors = errors,
        warnings = warnings,
        nv = nv,
        ne = ne,
        degree = degree,
        isolated = isolated,
        dangling = dangling,
        duplicates = duplicates,
        ncomponents = length(components),
        components = components,
        B1 = B1,
    )
end

# ---------- face-boundary helpers ----------

function _face_support_data(EV::Vector, face_edges::Vector{Int})
    edge_to_local = Dict{Int,Int}(e => k for (k,e) in enumerate(face_edges))
    incidence = Dict{Int, Vector{Tuple{Int,Int}}}()  # vertex => [(edge_id, other_vertex)]
    used_vertices_set = Set{Int}()

    for e in face_edges
        u, v = EV[e]
        push!(get!(incidence, u, Tuple{Int,Int}[]), (e, v))
        push!(get!(incidence, v, Tuple{Int,Int}[]), (e, u))
        push!(used_vertices_set, u)
        push!(used_vertices_set, v)
    end

    used_vertices = sort!(collect(used_vertices_set))
    return edge_to_local, incidence, used_vertices
end

function _connected_edge_components(EV::Vector, face_edges::Vector{Int}, incidence, used_vertices)
    vindex = Dict(v => i for (i,v) in enumerate(used_vertices))
    adjacency = [Int[] for _ in 1:length(used_vertices)]
    for v in used_vertices
        for (_, w) in incidence[v]
            push!(adjacency[vindex[v]], vindex[w])
        end
    end
    comps_local = graph_components(adjacency)
    comps_vertices = [used_vertices[idxs] for idxs in comps_local]

    comps_edges = Vector{Vector{Int}}(undef, length(comps_vertices))
    for (k, comp_vertices) in enumerate(comps_vertices)
        vertexset = Set(comp_vertices)
        comp_edges = Int[]
        for e in face_edges
            u, v = EV[e]
            if (u in vertexset) && (v in vertexset)
                push!(comp_edges, e)
            end
        end
        comps_edges[k] = comp_edges
    end
    return comps_vertices, comps_edges
end

"""
    _ordered_face_cycle(EV, face_edges)

Generalized exact combinatorial face-boundary extractor.
A face support may consist of one or more connected simple cycles (e.g. a planar
2-cell with holes). The routine returns coefficients aligned with `face_edges`.

If more than one connected component is present, the routine accepts the support
provided that every connected component is itself a simple cycle. Relative signs
between boundary components are assigned later by geometric nesting analysis.
"""
function _ordered_face_cycle(EV::Vector, face_edges::Vector{Int})
    m = length(face_edges)
    m >= 3 || return (ok=false, reason="face has fewer than 3 edges", coeffs=nothing,
                      used_vertices=Int[], ordered_edges=Vector{Vector{Int}}(),
                      component_reports=Any[])

    edge_to_local, incidence, used_vertices = _face_support_data(EV, face_edges)

    bad_degree_vertices = Int[]
    for v in used_vertices
        if length(get(incidence, v, Tuple{Int,Int}[])) != 2
            push!(bad_degree_vertices, v)
        end
    end
    if !isempty(bad_degree_vertices)
        return (ok=false,
                reason="face support is not a 2-regular subgraph; bad vertices: $(bad_degree_vertices)",
                coeffs=nothing,
                used_vertices=used_vertices,
                ordered_edges=Vector{Vector{Int}}(),
                component_reports=Any[])
    end

    comps_vertices, comps_edges = _connected_edge_components(EV, face_edges, incidence, used_vertices)

    coeffs = zeros(Int8, m)
    ordered_edges = Vector{Vector{Int}}()
    component_reports = Any[]

    for comp_edges in comps_edges
        comp_used_vertices = sort!(unique(vcat([collect(EV[e]) for e in comp_edges]...)))
        comp_incidence = Dict{Int, Vector{Tuple{Int,Int}}}()
        for e in comp_edges
            u, v = EV[e]
            push!(get!(comp_incidence, u, Tuple{Int,Int}[]), (e, v))
            push!(get!(comp_incidence, v, Tuple{Int,Int}[]), (e, u))
        end

        start = minimum(comp_used_vertices)
        nbrs = comp_incidence[start]
        first_edge, first_next = nbrs[1]

        comp_ordered_edges = Int[]
        traversed = Tuple{Int,Int,Int}[]  # (edge_id, from, to)
        visited_edges = Set{Int}()

        current_v = start
        next_v = first_next
        current_e = first_edge

        while true
            push!(comp_ordered_edges, current_e)
            push!(traversed, (current_e, current_v, next_v))
            push!(visited_edges, current_e)

            if next_v == start
                break
            end

            choices = comp_incidence[next_v]
            if choices[1][1] == current_e
                new_e, new_next = choices[2]
            else
                new_e, new_next = choices[1]
            end

            if new_e in visited_edges
                return (ok=false,
                        reason="premature edge repetition while traversing face boundary component",
                        coeffs=nothing,
                        used_vertices=used_vertices,
                        ordered_edges=ordered_edges,
                        component_reports=component_reports)
            end

            current_v, next_v, current_e = next_v, new_next, new_e
        end

        if length(visited_edges) != length(comp_edges)
            return (ok=false,
                    reason="face support contains extra edges outside a traversed cycle component",
                    coeffs=nothing,
                    used_vertices=used_vertices,
                    ordered_edges=ordered_edges,
                    component_reports=component_reports)
        end

        comp_vertices_ordered = [tr[2] for tr in traversed]
        comp_local_coeffs = Int8[]
        for (edge_id, fromv, tov) in traversed
            u, v = EV[edge_id]
            local_idx = edge_to_local[edge_id]
            s = (u == fromv && v == tov) ? Int8(+1) : Int8(-1)
            coeffs[local_idx] = s
            push!(comp_local_coeffs, s)
        end

        push!(ordered_edges, comp_ordered_edges)
        push!(component_reports, (
            ok = true,
            reason = "simple cycle component",
            vertices = comp_vertices_ordered,
            ordered_edges = comp_ordered_edges,
            local_coeffs = comp_local_coeffs,
        ))
    end

    reason = length(comps_edges) == 1 ? "single simple cycle" : "$(length(comps_edges)) simple boundary cycles"

    return (
        ok = true,
        reason = reason,
        coeffs = coeffs,
        used_vertices = used_vertices,
        ordered_edges = ordered_edges,
        component_reports = component_reports,
    )
end

function _find_face_plane_basis(V::AbstractMatrix, used_vertices::Vector{Int}; atol::Float64=1.0e-10)
    length(used_vertices) >= 3 || return nothing
    pts = [Vector{Float64}(V[:, v]) for v in used_vertices]
    p0 = pts[1]
    for i in 2:length(pts)-1
        a = pts[i] .- p0
        if norm(a) <= atol
            continue
        end
        for j in i+1:length(pts)
            b = pts[j] .- p0
            n = cross(a, b)
            if norm(n) > atol
                e1 = a / norm(a)
                n = n / norm(n)
                e2 = cross(n, e1)
                return (p0=p0, e1=e1, e2=e2, n=n)
            end
        end
    end
    return nothing
end

_project_2d(p::AbstractVector, p0, e1, e2) = [dot(p .- p0, e1), dot(p .- p0, e2)]

function _polygon_signed_area(poly::Vector{Vector{Float64}})
    n = length(poly)
    A = 0.0
    for i in 1:n
        p = poly[i]
        q = poly[i == n ? 1 : i+1]
        A += p[1]*q[2] - q[1]*p[2]
    end
    return 0.5*A
end

function _point_in_polygon(pt::Vector{Float64}, poly::Vector{Vector{Float64}}; atol::Float64=1.0e-12)
    x, y = pt
    inside = false
    n = length(poly)
    for i in 1:n
        p = poly[i]
        q = poly[i == n ? 1 : i+1]
        xi, yi = p
        xj, yj = q

        # Boundary check (collinearity + bounding box)
        dx1, dy1 = x - xi, y - yi
        dx2, dy2 = xj - xi, yj - yi
        crossv = dx1*dy2 - dy1*dx2
        if abs(crossv) <= atol && min(xi,xj)-atol <= x <= max(xi,xj)+atol && min(yi,yj)-atol <= y <= max(yi,yj)+atol
            return true
        end

        intersect = ((yi > y) != (yj > y)) &&
                    (x < (xj - xi) * (y - yi) / (yj - yi + eps()) + xi)
        if intersect
            inside = !inside
        end
    end
    return inside
end

function _sample_interior_point(poly::Vector{Vector{Float64}}, area::Float64)
    n = length(poly)
    if n == 0
        return [0.0, 0.0]
    elseif n == 1
        return copy(poly[1])
    end
    p1 = poly[1]
    p2 = poly[2]
    t = p2 .- p1
    scale = max(maximum(abs, reduce(vcat, poly)), 1.0)
    ε = 1.0e-6 * scale
    left_normal = [-t[2], t[1]]
    norm(left_normal) <= eps() && return 0.5 .* (p1 .+ p2)
    inward = area >= 0 ? left_normal : -left_normal
    inward ./= norm(inward)
    return 0.5 .* (p1 .+ p2) .+ ε .* inward
end

function _assign_component_nesting_signs(V::AbstractMatrix, used_vertices::Vector{Int}, component_reports; atol::Float64=1.0e-10)
    basis = _find_face_plane_basis(V, used_vertices; atol=atol)
    basis === nothing && return (ok=false, reason="cannot determine a supporting plane for face boundary components", multipliers=nothing)

    p0, e1, e2 = basis.p0, basis.e1, basis.e2
    ncomp = length(component_reports)
    polygons = Vector{Vector{Vector{Float64}}}(undef, ncomp)
    areas = zeros(Float64, ncomp)
    sample_points = Vector{Vector{Float64}}(undef, ncomp)

    for k in 1:ncomp
        verts = component_reports[k].vertices
        poly = [_project_2d(Vector{Float64}(V[:, v]), p0, e1, e2) for v in verts]
        A = _polygon_signed_area(poly)
        abs(A) > atol || return (ok=false, reason="degenerate projected boundary component with near-zero area", multipliers=nothing)
        polygons[k] = poly
        areas[k] = A
        sample_points[k] = _sample_interior_point(poly, A)
    end

    depths = zeros(Int, ncomp)
    for i in 1:ncomp
        for j in 1:ncomp
            i == j && continue
            if _point_in_polygon(sample_points[i], polygons[j])
                depths[i] += 1
            end
        end
    end

    roots = findall(==(0), depths)
    length(roots) == 1 || return (ok=false, reason="boundary components do not define one connected planar face region (root components: $(roots))", multipliers=nothing)

    multipliers = ones(Int8, ncomp)
    for i in 1:ncomp
        desired_sign = iseven(depths[i]) ? 1.0 : -1.0
        multipliers[i] = sign(areas[i]) == sign(desired_sign) ? Int8(+1) : Int8(-1)
    end

    return (ok=true, reason="success", multipliers=multipliers, depths=depths, areas=areas)
end

"""
    test_face_cycle(V, EV, face_edges)

Exact combinatorial/geometric validator for one face support.
Accepts one or more disjoint simple cycle components, provided they define the
boundary of one connected planar face region (outer loop plus holes).
"""
function test_face_cycle(V::AbstractMatrix, EV::Vector, face_edges::Vector{Int})
    result = _ordered_face_cycle(EV, face_edges)
    if !result.ok
        return (
            ok = false,
            reason = result.reason,
            coeffs = nothing,
            used_vertices = result.used_vertices,
            ordered_edges = result.ordered_edges,
            ncomponents = length(result.component_reports),
            component_reports = result.component_reports,
        )
    end

    coeffs = copy(result.coeffs)
    ncomp = length(result.component_reports)
    nesting = nothing

    if ncomp > 1
        nesting = _assign_component_nesting_signs(V, result.used_vertices, result.component_reports)
        if !nesting.ok
            return (
                ok = false,
                reason = nesting.reason,
                coeffs = nothing,
                used_vertices = result.used_vertices,
                ordered_edges = result.ordered_edges,
                ncomponents = ncomp,
                component_reports = result.component_reports,
            )
        end
        for (k, comp) in enumerate(result.component_reports)
            mult = nesting.multipliers[k]
            for (j, edge_id) in enumerate(comp.ordered_edges)
                local_idx = findfirst(==(edge_id), face_edges)
                coeffs[local_idx] = Int8(mult * comp.local_coeffs[j])
            end
        end
    end

    return (
        ok = true,
        reason = ncomp == 1 ? "single simple cycle" : "$(ncomp) boundary cycles (face with holes)",
        coeffs = coeffs,
        used_vertices = result.used_vertices,
        ordered_edges = result.ordered_edges,
        ncomponents = ncomp,
        component_reports = result.component_reports,
        nesting = nesting,
    )
end

"""
    validate_face_cycles(V, EV, FE)

Check every face support in FE and, if all pass, build a provisional EF0 matrix.
Each face column is individually a closed 1-cycle, and faces with holes are
accepted when their support is the disjoint union of simple boundary cycles.
"""
function validate_face_cycles(V::AbstractMatrix, EV::Vector, FE::Vector)
    nf = length(FE)
    ne = length(EV)

    reports = Vector{Any}(undef, nf)
    I = Int[]
    J = Int[]
    X = Int8[]

    ok = true
    for f in 1:nf
        report = test_face_cycle(V, EV, FE[f])
        reports[f] = report
        if !report.ok
            ok = false
            continue
        end
        for (k,e) in enumerate(FE[f])
            push!(I, e)
            push!(J, f)
            push!(X, report.coeffs[k])
        end
    end

    EF0 = ok ? sparse(I, J, X, ne, nf) : nothing
    failed_faces = findall(r -> !r.ok, reports)
    return (ok=ok, reports=reports, failed_faces=failed_faces, EF0=EF0)
end

"""
    edge_face_multiplicity(EF)

Number of incident faces for each edge row.
"""
function edge_face_multiplicity(EF::SparseMatrixCSC)
    ne, _ = size(EF)
    mult = zeros(Int, ne)
    for e in 1:ne
        _, vals = findnz(EF[e, :])
        mult[e] = length(vals)
    end
    return mult
end

"""
    _row_data(EF0)

Extract row supports and coefficients from EF0.
Returns vectors `rowfaces`, `rowvals`, one entry per edge row.
"""
function _row_data(EF0::SparseMatrixCSC)
    ne, _ = size(EF0)
    rowfaces = Vector{Vector{Int}}(undef, ne)
    rowvals  = Vector{Vector{Int8}}(undef, ne)
    for e in 1:ne
        Js, Vs = findnz(EF0[e, :])
        rowfaces[e] = collect(Js)
        rowvals[e] = collect(Vs)
    end
    return rowfaces, rowvals
end

"""
    _balance_possible(sum_now, k_remaining)

Pruning predicate for row-balance equations with remaining `k_remaining` variables,
where each remaining term contributes ±1.
Target is final sum zero.
"""
function _balance_possible(sum_now::Int, k_remaining::Int)
    abs(sum_now) <= k_remaining || return false
    return iseven(k_remaining - abs(sum_now))
end

"""
    solve_face_flip_balance(EF0)

Given a provisional EF0 whose columns are individually closed cycles, try to find
face signs s_f ∈ {±1} such that every edge row sums to zero:
    Σ_f EF0[e,f] * s_f = 0

This accepts any even edge multiplicity. Hence it validates algebraic closedness
of the 2-skeleton even in non-manifold cases.
"""
function solve_face_flip_balance(EF0::SparseMatrixCSC)
    ne, nf = size(EF0)
    rowfaces, rowvals = _row_data(EF0)
    multiplicity = [length(rowfaces[e]) for e in 1:ne]
    odd_rows = findall(x -> isodd(x), multiplicity)

    if !isempty(odd_rows)
        return (
            ok = false,
            reason = "rows with odd edge-face multiplicity: $(odd_rows)",
            face_signs = zeros(Int8, nf),
            EF = nothing,
            multiplicity = multiplicity,
            odd_rows = odd_rows,
            nonzero_row_sums = Int[],
            explored_nodes = 0,
        )
    end

    facerows = [Int[] for _ in 1:nf]
    faceslots = [Int[] for _ in 1:nf]
    for e in 1:ne
        for (slot, f) in enumerate(rowfaces[e])
            push!(facerows[f], e)
            push!(faceslots[f], slot)
        end
    end

    assigned = fill(Int8(0), nf)
    rowsum = zeros(Int, ne)
    rowunassigned = copy(multiplicity)
    explored = Ref(0)

    function choose_face()
        bestf = 0
        bestscore = -1
        for f in 1:nf
            if assigned[f] != 0
                continue
            end
            score = 0
            for e in facerows[f]
                score += 10 - min(rowunassigned[e], 9)
                score += abs(rowsum[e])
            end
            if score > bestscore
                bestscore = score
                bestf = f
            end
        end
        return bestf
    end

    function assign_face!(f::Int, s::Int8)
        assigned[f] = s
        for (idx, e) in enumerate(facerows[f])
            slot = faceslots[f][idx]
            rowsum[e] += Int(rowvals[e][slot]) * Int(s)
            rowunassigned[e] -= 1
        end
    end

    function unassign_face!(f::Int, s::Int8)
        for (idx, e) in enumerate(facerows[f])
            slot = faceslots[f][idx]
            rowsum[e] -= Int(rowvals[e][slot]) * Int(s)
            rowunassigned[e] += 1
        end
        assigned[f] = 0
    end

    function feasible_rows_after(f::Int)
        for e in facerows[f]
            _balance_possible(rowsum[e], rowunassigned[e]) || return false
        end
        return true
    end

    function all_rows_zero()
        for e in 1:ne
            rowsum[e] == 0 || return false
        end
        return true
    end

    function search!()
        explored[] += 1
        f = choose_face()
        if f == 0
            return all_rows_zero()
        end

        for s in (Int8(+1), Int8(-1))
            assign_face!(f, s)
            if feasible_rows_after(f)
                globally_ok = true
                for e in 1:ne
                    if rowunassigned[e] == 0
                        if rowsum[e] != 0
                            globally_ok = false
                            break
                        end
                    elseif !_balance_possible(rowsum[e], rowunassigned[e])
                        globally_ok = false
                        break
                    end
                end
                if globally_ok && search!()
                    return true
                end
            end
            unassign_face!(f, s)
        end
        return false
    end

    if nf > 0
        assign_face!(1, Int8(+1))
        root_ok = true
        for e in 1:ne
            if !_balance_possible(rowsum[e], rowunassigned[e])
                root_ok = false
                break
            end
        end
        solved = root_ok && search!()
        if !solved
            unassign_face!(1, Int8(+1))
        end
    else
        solved = true
    end

    if !solved
        return (
            ok = false,
            reason = "no global face-sign assignment makes every edge row sum to zero",
            face_signs = copy(assigned),
            EF = nothing,
            multiplicity = multiplicity,
            odd_rows = odd_rows,
            nonzero_row_sums = findall(!=(0), rowsum),
            explored_nodes = explored[],
        )
    end

    I, J, X = findnz(EF0)
    Xnew = Int8[]
    for k in eachindex(X)
        push!(Xnew, Int8(X[k] * assigned[J[k]]))
    end
    EF = sparse(I, J, Xnew, ne, nf)

    final_rowsums = [sum(last(findnz(EF[e, :]))) for e in 1:ne]
    bad_rows = findall(!=(0), final_rowsums)
    reason = isempty(bad_rows) ? "success" : "internal error: nonzero row sums remain"

    return (
        ok = isempty(bad_rows),
        reason = reason,
        face_signs = assigned,
        EF = isempty(bad_rows) ? EF : nothing,
        multiplicity = multiplicity,
        odd_rows = odd_rows,
        nonzero_row_sums = bad_rows,
        explored_nodes = explored[],
    )
end

"""
    validate_2skeleton(V, EV, FE; manifold_check=false, expect_no_isolated=true)

Main validation routine for the global embedded 2-skeleton.

Pipeline:
1. validate EV and build B1;
2. validate each face support as one or more exact boundary cycles;
3. build a provisional EF0 (each face closed independently);
4. solve the global face-flip problem so that every edge row of EF sums to zero;
5. optionally report manifoldness diagnostics (edge multiplicity == 2).
"""
function validate_2skeleton(V::AbstractMatrix, EV::Vector, FE::Vector;
                            manifold_check::Bool=false,
                            expect_no_isolated::Bool=true)
    evdiag = validate_EV(V, EV; expect_no_isolated=expect_no_isolated)
    B1 = evdiag.B1
    if !evdiag.ok
        return (
            ok = false,
            stage = :EV,
            EV = evdiag,
            ev = evdiag,
            facecycles = nothing,
            faces = nothing,
            balance = nothing,
            orientation = nothing,
            manifold = nothing,
            B1 = nothing,
            EF0 = nothing,
            EF = nothing,
        )
    end

    facediag = validate_face_cycles(V, EV, FE)
    if !facediag.ok
        return (
            ok = false,
            stage = :facecycles,
            EV = evdiag,
            ev = evdiag,
            facecycles = facediag,
            faces = facediag,
            balance = nothing,
            orientation = nothing,
            manifold = nothing,
            B1 = B1,
            EF0 = nothing,
            EF = nothing,
        )
    end

    balancediag = solve_face_flip_balance(facediag.EF0)
    if !balancediag.ok
        manifold = (
            ok = false,
            multiplicity = balancediag.multiplicity,
            nonmanifold_edges = findall(!=(2), balancediag.multiplicity),
        )
        return (
            ok = false,
            stage = :balance,
            EV = evdiag,
            ev = evdiag,
            facecycles = facediag,
            faces = facediag,
            balance = balancediag,
            orientation = balancediag,
            manifold = manifold,
            B1 = B1,
            EF0 = facediag.EF0,
            EF = nothing,
        )
    end

    EF = balancediag.EF
    chain_mat = B1 * EF
    chain_ok = maximum(abs, Matrix(chain_mat)) == 0
    row_sums = vec(sum(EF, dims=2))
    balanced_rows_ok = all(row_sums .== 0)

    mult = balancediag.multiplicity
    nonmanifold_edges = findall(!=(2), mult)
    manifold = (
        ok = isempty(nonmanifold_edges),
        multiplicity = mult,
        nonmanifold_edges = nonmanifold_edges,
    )

    algebraic_ok = evdiag.ok && facediag.ok && balancediag.ok && chain_ok && balanced_rows_ok
    final_ok = manifold_check ? (algebraic_ok && manifold.ok) : algebraic_ok

    return (
        ok = final_ok,
        stage = :done,
        EV = evdiag,
        ev = evdiag,
        facecycles = facediag,
        faces = facediag,
        balance = balancediag,
        orientation = balancediag,
        manifold = manifold,
        B1 = B1,
        EF0 = facediag.EF0,
        EF = EF,
        chain_ok = chain_ok,
        balanced_rows_ok = balanced_rows_ok,
        row_sums = row_sums,
    )
end

end # module

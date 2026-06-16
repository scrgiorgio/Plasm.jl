module BlockPairCongruence

using SparseArrays
using LinearAlgebra

export block_pairs_from_FEorig,
       edge_incident_blocks,
       edge_segment_from_lar,
       boundary_vertices_on_original_edge,
       boundary_edges_on_original_edge,
       compare_edge_subdivision,
       validate_all_original_edges,
       pair_congruence_report

export repair_all_original_edges,
       assemble_repaired_blocks,
       repair_report

# ------------------------------------------------------------
# Basic sparse helpers
# ------------------------------------------------------------

"""
    _col_support(M, j)

Return row indices of the nonzeros in column `j` of a sparse matrix.
"""
function _col_support(M::SparseMatrixCSC, j::Int)
    rows, _ = findnz(M[:, j])
    return collect(rows)
end

"""
    _row_support(M, i)

Return column indices of the nonzeros in row `i` of a sparse matrix.
Works directly on SparseMatrixCSC.
"""
function _row_support(M::SparseMatrixCSC, i::Int)
    cols = Int[]
    @inbounds for j in 1:size(M, 2)
        for p in nzrange(M, j)
            if M.rowval[p] == i
                push!(cols, j)
            end
        end
    end
    return cols
end

# ------------------------------------------------------------
# Geometry helpers
# ------------------------------------------------------------

"""
    _point_col(W, i)

Return point `i` as a dense vector from a d×n matrix W.
"""
_point_col(W::AbstractMatrix, i::Int) = Vector(W[:, i])

"""
    _segment_parameter(x, p, q)

Return the affine parameter t of point x projected onto segment pq:
    x ≈ p + t (q-p)
"""
function _segment_parameter(x::AbstractVector, p::AbstractVector, q::AbstractVector)
    d = q - p
    den = dot(d, d)
    den == 0.0 && error("degenerate original edge segment")
    return dot(x - p, d) / den
end

"""
    _point_line_distance(x, p, q)

Euclidean distance from x to the line through p,q.
"""
function _point_line_distance(x::AbstractVector, p::AbstractVector, q::AbstractVector)
    d = q - p
    den = norm(d)
    den == 0.0 && return norm(x - p)
    t = dot(x - p, d) / dot(d, d)
    proj = p + t * d
    return norm(x - proj)
end

"""
    _point_on_segment(x, p, q; atol)

Test if point x lies on the closed segment pq within tolerance.
Returns `(ok, t, dist)`.
"""
function _point_on_segment(x::AbstractVector, p::AbstractVector, q::AbstractVector; atol::Float64=1e-7)
    t = _segment_parameter(x, p, q)
    dist = _point_line_distance(x, p, q)
    ok = (-atol <= t <= 1 + atol) && (dist <= atol)
    return (ok=ok, t=t, dist=dist)
end

"""
    _normalize_t(t; atol)

Clamp near-endpoint values to 0 or 1.
"""
function _normalize_t(t::Real; atol::Float64=1e-7)
    if abs(t) <= atol
        return 0.0
    elseif abs(t - 1.0) <= atol
        return 1.0
    else
        return Float64(t)
    end
end

"""
    _unique_sorted(vals; atol)

Sort numeric values and merge values within tolerance.
"""
function _unique_sorted(vals::Vector{Float64}; atol::Float64=1e-7)
    isempty(vals) && return Float64[]
    v = sort(vals)
    out = Float64[v[1]]
    for x in v[2:end]
        if abs(x - out[end]) > atol
            push!(out, x)
        end
    end
    return out
end

"""
    _intervals_from_tvals(tvals; atol)

Given sorted unique parameter values, return consecutive intervals.
"""
function _intervals_from_tvals(tvals::Vector{Float64}; atol::Float64=1e-7)
    if length(tvals) < 2
        return Tuple{Float64,Float64}[]
    end
    out = Tuple{Float64,Float64}[]
    for i in 1:length(tvals)-1
        a, b = tvals[i], tvals[i+1]
        if b - a > atol
            push!(out, (a, b))
        end
    end
    return out
end

"""
    _same_tsets(a, b; atol)

Compare two sorted parameter sets up to tolerance.
"""
function _same_tsets(a::Vector{Float64}, b::Vector{Float64}; atol::Float64=1e-7)
    length(a) == length(b) || return false
    all(abs(a[i] - b[i]) <= atol for i in eachindex(a))
end

# ------------------------------------------------------------
# Original topology from FE_orig and lar
# ------------------------------------------------------------

"""
    edge_incident_blocks(FE_orig, e)

Return the two original face/block indices incident on original edge e.
"""
function edge_incident_blocks(FE_orig::SparseMatrixCSC, e::Int)
    rows = sort(_col_support(FE_orig, e))
    return rows
end

"""
    block_pairs_from_FEorig(FE_orig)

Return a vector of named tuples, one per original edge, containing:
- original edge id
- the pair of incident block ids
"""
function block_pairs_from_FEorig(FE_orig::SparseMatrixCSC)
    pairs = NamedTuple[]
    for e in 1:size(FE_orig, 2)
        blocks = edge_incident_blocks(FE_orig, e)
        push!(pairs, (
            edge = e,
            blocks = blocks,
            ok = (length(blocks) == 2),
        ))
    end
    return pairs
end

"""
    edge_segment_from_lar(lar, e)

Return endpoints and coordinates of original edge e from lar.
Assumes:
- `lar.V` is d×n
- `lar.C[:EV]` is the original edge-vertex list
"""
function edge_segment_from_lar(lar, e::Int)
    EVorig = lar.C[:EV]
    u, v = EVorig[e]
    p = Vector(lar.V[:, u])
    q = Vector(lar.V[:, v])
    return (u=u, v=v, p=p, q=q)
end

# ------------------------------------------------------------
# Local block boundary support restricted to one original edge
# ------------------------------------------------------------

"""
    _local_EV(block)

Build local edge-vertex list from block.Delta0 (ne × nv).
"""
function _local_EV(block)
    Delta0 = block.Delta0
    EV = Vector{Vector{Int}}(undef, size(Delta0, 1))
    for e in 1:size(Delta0, 1)
        EV[e] = sort(_row_support(Delta0, e))
    end
    return EV
end

"""
    _boundary_edge_ids(block)

Local edge ids used by exactly one local arranged face.
"""
function _boundary_edge_ids(block)
    Delta1 = block.Delta1
    ne = size(Delta1, 2)
    b = Int[]
    for e in 1:ne
        rows, _ = findnz(Delta1[:, e])
        if length(rows) == 1
            push!(b, e)
        end
    end
    return b
end

"""
    boundary_vertices_on_original_edge(block, p, q; atol)

Return all local boundary vertices of `block` lying on original segment pq.
Also returns their affine parameters.
"""
function boundary_vertices_on_original_edge(block, p::AbstractVector, q::AbstractVector; atol::Float64=1e-7)
    W = block.W
    EV = _local_EV(block)
    bedges = _boundary_edge_ids(block)

    bverts = sort(unique(v for e in bedges for v in EV[e]))

    # allow temporary/repaired blocks that do not carry global ids
    gv_ids = hasproperty(block, :global_vertex_ids) ? collect(block.global_vertex_ids) : collect(1:size(W, 2))

    local_vertices = Int[]
    global_vertices = Int[]
    tvals = Float64[]
    coords = Vector{Float64}[]

    for v in bverts
        x = _point_col(W, v)
        test = _point_on_segment(x, p, q; atol=atol)
        if test.ok
            push!(local_vertices, v)
            push!(global_vertices, gv_ids[v])
            push!(tvals, _normalize_t(test.t; atol=atol))
            push!(coords, x)
        end
    end

    ord = sortperm(tvals)
    return (
        local_vertices = local_vertices[ord],
        global_vertices = global_vertices[ord],
        tvals = tvals[ord],
        coords = coords[ord],
    )
end

"""
    boundary_edges_on_original_edge(block, p, q; atol)

Return all local boundary edges of `block` whose both endpoints lie on original segment pq.
Also returns parameter intervals.
"""
function boundary_edges_on_original_edge(block, p::AbstractVector, q::AbstractVector; atol::Float64=1e-7)
    W = block.W
    EV = _local_EV(block)
    bedges = _boundary_edge_ids(block)

    # allow temporary/repaired blocks that do not carry global ids
    gv_ids = hasproperty(block, :global_vertex_ids) ? collect(block.global_vertex_ids) : collect(1:size(W, 2))
    ge_ids = hasproperty(block, :global_edge_ids) ? collect(block.global_edge_ids) : collect(1:size(block.Delta0, 1))

    local_edges = Int[]
    global_edges = Int[]
    intervals = Tuple{Float64,Float64}[]
    local_endpoints = Vector{Tuple{Int,Int}}()
    global_endpoints = Vector{Tuple{Int,Int}}()

    for e in bedges
        a, b = EV[e]
        xa = _point_col(W, a)
        xb = _point_col(W, b)

        ta = _point_on_segment(xa, p, q; atol=atol)
        tb = _point_on_segment(xb, p, q; atol=atol)

        if ta.ok && tb.ok
            t1 = _normalize_t(ta.t; atol=atol)
            t2 = _normalize_t(tb.t; atol=atol)
            if abs(t2 - t1) > atol
                lo, hi = minmax(t1, t2)
                push!(local_edges, e)
                push!(global_edges, ge_ids[e])
                push!(intervals, (lo, hi))
                push!(local_endpoints, (a, b))
                push!(global_endpoints, (gv_ids[a], gv_ids[b]))
            end
        end
    end

    ord = sortperm(first.(intervals))
    return (
        local_edges = local_edges[ord],
        global_edges = global_edges[ord],
        intervals = intervals[ord],
        local_endpoints = local_endpoints[ord],
        global_endpoints = global_endpoints[ord],
    )
end
# ------------------------------------------------------------
# Pairwise comparison along one original edge
# ------------------------------------------------------------

"""
    compare_edge_subdivision(block1, block2, p, q; atol)

Compare the induced boundary subdivision of original edge segment pq
between two adjacent blocks.
"""
function compare_edge_subdivision(block1, block2, p::AbstractVector, q::AbstractVector; atol::Float64=1e-7)
    V1 = boundary_vertices_on_original_edge(block1, p, q; atol=atol)
    V2 = boundary_vertices_on_original_edge(block2, p, q; atol=atol)

    E1 = boundary_edges_on_original_edge(block1, p, q; atol=atol)
    E2 = boundary_edges_on_original_edge(block2, p, q; atol=atol)

    t1 = _unique_sorted(copy(V1.tvals); atol=atol)
    t2 = _unique_sorted(copy(V2.tvals); atol=atol)

    ints1 = _intervals_from_tvals(t1; atol=atol)
    ints2 = _intervals_from_tvals(t2; atol=atol)

    same_vertices = _same_tsets(t1, t2; atol=atol)
    same_intervals = (length(ints1) == length(ints2)) &&
                     all(abs(ints1[i][1] - ints2[i][1]) <= atol &&
                         abs(ints1[i][2] - ints2[i][2]) <= atol
                         for i in eachindex(ints1))

    only1 = Float64[]
    only2 = Float64[]
    for x in t1
        if !any(abs(x - y) <= atol for y in t2)
            push!(only1, x)
        end
    end
    for y in t2
        if !any(abs(y - x) <= atol for x in t1)
            push!(only2, y)
        end
    end

    ok = same_vertices && same_intervals

    return (
        ok = ok,
        reason = ok ? "matching subdivision along original edge" :
                      "mismatching boundary subdivision along original edge",
        block1_vertices = V1,
        block2_vertices = V2,
        block1_edges = E1,
        block2_edges = E2,
        tvals1 = t1,
        tvals2 = t2,
        intervals1 = ints1,
        intervals2 = ints2,
        only_in_block1 = only1,
        only_in_block2 = only2,
    )
end

# ------------------------------------------------------------
# Global driver over original edges
# ------------------------------------------------------------

"""
    validate_all_original_edges(blocks, FE_orig, lar; atol)

For each original edge of the input complex:
- read its incident original face/block pair from FE_orig,
- compare the two local block boundary subdivisions along that edge.

Returns a global report.
"""
function validate_all_original_edges(blocks::Vector,
                                     FE_orig::SparseMatrixCSC,
                                     lar;
                                     atol::Float64=1e-7)

    pairs = block_pairs_from_FEorig(FE_orig)
    reports = Any[]
    failed_edges = Int[]
    malformed_edges = Int[]

    for item in pairs
        e = item.edge
        incident = item.blocks

        if length(incident) != 2
            push!(malformed_edges, e)
            push!(reports, (
                ok = false,
                edge = e,
                reason = "original edge does not have exactly two incident blocks",
                incident_blocks = incident,
            ))
            continue
        end

        σ1, σ2 = incident
        block1 = blocks[σ1]
        block2 = blocks[σ2]
        seg = edge_segment_from_lar(lar, e)

        cmp = compare_edge_subdivision(block1, block2, seg.p, seg.q; atol=atol)

        rep = (
            ok = cmp.ok,
            edge = e,
            incident_blocks = (σ1, σ2),
            original_edge_vertices = (seg.u, seg.v),
            p = seg.p,
            q = seg.q,
            comparison = cmp,
        )

        if !cmp.ok
            push!(failed_edges, e)
        end
        push!(reports, rep)
    end

    return (
        ok = isempty(failed_edges) && isempty(malformed_edges),
        nedges = size(FE_orig, 2),
        failed_edges = failed_edges,
        malformed_edges = malformed_edges,
        reports = reports,
    )
end

# ------------------------------------------------------------
# Text report
# ------------------------------------------------------------

"""
    pair_congruence_report(summary)

Compact textual report for pairwise boundary congruence along original edges.
"""
function pair_congruence_report(summary)
    lines = String[]
    push!(lines, "Block pair congruence summary")
    push!(lines, "  number of original edges checked : $(summary.nedges)")
    push!(lines, "  all pair checks passed           : $(summary.ok)")
    push!(lines, "  malformed original edges         : $(summary.malformed_edges)")
    push!(lines, "  failed original edges            : $(summary.failed_edges)")

    for rep in summary.reports
        if get(rep, :ok, false)
            continue
        end

        if !haskey(rep, :comparison)
            push!(lines, "")
            push!(lines, "  edge $(rep.edge): ok = false")
            push!(lines, "    reason = $(rep.reason)")
            push!(lines, "    incident blocks = $(rep.incident_blocks)")
            continue
        end

        cmp = rep.comparison
        push!(lines, "")
        push!(lines, "  edge $(rep.edge): ok = false")
        push!(lines, "    original edge vertices = $(rep.original_edge_vertices)")
        push!(lines, "    incident blocks        = $(rep.incident_blocks)")
        push!(lines, "    reason                 = $(cmp.reason)")
        push!(lines, "    block1 t-values        = $(cmp.tvals1)")
        push!(lines, "    block2 t-values        = $(cmp.tvals2)")
        push!(lines, "    only in block1         = $(cmp.only_in_block1)")
        push!(lines, "    only in block2         = $(cmp.only_in_block2)")
        push!(lines, "    block1 intervals       = $(cmp.intervals1)")
        push!(lines, "    block2 intervals       = $(cmp.intervals2)")
    end

    return join(lines, "\n")
end

# ------------------------------------------------------------
# REPAIR SUPPORT
# ------------------------------------------------------------

"""
    _find_existing_vertex(W, x; atol)

Return the local vertex id whose coordinates match `x` within tolerance,
or 0 if none exists.
"""
function _find_existing_vertex(W::AbstractMatrix, x::AbstractVector; atol::Float64=1e-7)
    for v in 1:size(W, 2)
        if norm(W[:, v] - x) <= atol
            return v
        end
    end
    return 0
end

"""
    _append_vertex(W, x)

Append point x as a new column to W.
"""
function _append_vertex(W::Matrix{Float64}, x::AbstractVector)
    return hcat(W, reshape(Vector{Float64}(x), :, 1))
end

"""
    _edge_rows(Delta0)

Extract each row of Delta0 as a pair of (vertex, sign) tuples.
Assumes each edge row has exactly two nonzeros.
"""
function _edge_rows(Delta0::SparseMatrixCSC{Int8, Int64})
    rows = Vector{Vector{Tuple{Int,Int8}}}(undef, size(Delta0, 1))
    for e in 1:size(Delta0, 1)
        cols = _row_support(Delta0, e)
        vals = Int8[]
        for c in cols
            push!(vals, Delta0[e, c])
        end
        rows[e] = [(cols[i], vals[i]) for i in eachindex(cols)]
    end
    return rows
end

"""
    _face_rows(Delta1)

Extract each row of Delta1 as a vector of (edge, sign) tuples.
"""
function _face_rows(Delta1::SparseMatrixCSC{Int8, Int64})
    rows = Vector{Vector{Tuple{Int,Int8}}}(undef, size(Delta1, 1))
    for f in 1:size(Delta1, 1)
        edges = _row_support(Delta1, f)
        vals = Int8[]
        for e in edges
            push!(vals, Delta1[f, e])
        end
        rows[f] = [(edges[i], vals[i]) for i in eachindex(edges)]
    end
    return rows
end

"""
    _build_Delta0(edge_rows, nv)

Rebuild sparse Delta0 from row representation.
"""
function _build_Delta0(edge_rows::Vector{Vector{Tuple{Int,Int8}}}, nv::Int)
    I = Int[]
    J = Int[]
    V = Int8[]
    for e in 1:length(edge_rows)
        for (v, s) in edge_rows[e]
            push!(I, e)
            push!(J, v)
            push!(V, s)
        end
    end
    return sparse(I, J, V, length(edge_rows), nv)
end

"""
    _build_Delta1(face_rows, ne)

Rebuild sparse Delta1 from row representation.
"""
function _build_Delta1(face_rows::Vector{Vector{Tuple{Int,Int8}}}, ne::Int)
    I = Int[]
    J = Int[]
    V = Int8[]
    for f in 1:length(face_rows)
        for (e, s) in face_rows[f]
            push!(I, f)
            push!(J, e)
            push!(V, s)
        end
    end
    return sparse(I, J, V, length(face_rows), ne)
end

"""
    _tail_head(edge_row)

Given one Delta0 row as [(v1,s1),(v2,s2)], return (tail, head).
Uses -1 = tail, +1 = head when available.
"""
function _tail_head(edge_row::Vector{Tuple{Int,Int8}})
    @assert length(edge_row) == 2
    (v1, s1), (v2, s2) = edge_row
    if s1 == -1 && s2 == 1
        return (v1, v2)
    elseif s1 == 1 && s2 == -1
        return (v2, v1)
    else
        # fallback if signs are malformed
        return (v1, v2)
    end
end

"""
    _replace_edge_indices(face_rows, e, k)

In every face row:
- replace old edge index e by k consecutive edges e:e+k-1 with same sign,
- shift all later edge indices by k-1.
"""
function _replace_edge_indices(face_rows::Vector{Vector{Tuple{Int,Int8}}}, e::Int, k::Int)
    out = Vector{Vector{Tuple{Int,Int8}}}(undef, length(face_rows))
    shift = k - 1
    for f in 1:length(face_rows)
        newrow = Tuple{Int,Int8}[]
        for (edgeid, s) in face_rows[f]
            if edgeid < e
                push!(newrow, (edgeid, s))
            elseif edgeid == e
                for ee in e:(e + k - 1)
                    push!(newrow, (ee, s))
                end
            else
                push!(newrow, (edgeid + shift, s))
            end
        end
        out[f] = newrow
    end
    return out
end

"""
    _repair_block_edge_to_tvals(block, p, q, target_tvals; atol)

Refine one block so that its boundary subdivision along original edge pq
contains all target_tvals.
"""
function _repair_block_edge_to_tvals(block, p::AbstractVector, q::AbstractVector,
                                     target_tvals::Vector{Float64}; atol::Float64=1e-7)

    W = Matrix{Float64}(block.W)
    Delta0 = SparseMatrixCSC{Int8,Int64}(block.Delta0)
    Delta1 = SparseMatrixCSC{Int8,Int64}(block.Delta1)

    changed = false
    progress = true
    maxiter = 100

    iter = 0
    while progress && iter < maxiter
        iter += 1
        progress = false

        Vrep = boundary_vertices_on_original_edge((W=W, Delta0=Delta0, Delta1=Delta1), p, q; atol=atol)
        Erep = boundary_edges_on_original_edge((W=W, Delta0=Delta0, Delta1=Delta1), p, q; atol=atol)

        have = _unique_sorted(copy(Vrep.tvals); atol=atol)
        missing = Float64[t for t in target_tvals if !any(abs(t - h) <= atol for h in have)]

        isempty(missing) && break

        # Try to split one boundary edge containing at least one missing t
        for idx in eachindex(Erep.local_edges)
            e = Erep.local_edges[idx]
            (a, b) = Erep.intervals[idx]
            mids = sort([t for t in missing if a + atol < t < b - atol])
            isempty(mids) && continue

            # Current combinatorial rows
            erows = _edge_rows(Delta0)
            frows = _face_rows(Delta1)

            oldrow = erows[e]
            tail, head = _tail_head(oldrow)

            xtail = _point_col(W, tail)
            xhead = _point_col(W, head)

            ttail = _normalize_t(_segment_parameter(xtail, p, q); atol=atol)
            thead = _normalize_t(_segment_parameter(xhead, p, q); atol=atol)

            # create/get interior vertices
            mids_sorted = sort(mids)
            mids_vids = Int[]
            mids_geom = Vector{Float64}[]
            for t in mids_sorted
                x = p + t * (q - p)
                vid = _find_existing_vertex(W, x; atol=atol)
                if vid == 0
                    W = _append_vertex(W, x)
                    vid = size(W, 2)
                end
                push!(mids_vids, vid)
                push!(mids_geom, x)
            end

            # order vertices along edge orientation tail -> head
            if ttail < thead
                chain = [tail; mids_vids; head]
            else
                chain = [tail; reverse(mids_vids); head]
            end

            # replace old edge row by chain of new oriented subedges
            new_edge_rows = copy(erows)
            segrows = Vector{Vector{Tuple{Int,Int8}}}()
            for i in 1:length(chain)-1
                push!(segrows, [(chain[i], Int8(-1)), (chain[i+1], Int8(1))])
            end

            # replace row e and insert the others after e
            new_edge_rows[e] = segrows[1]
            if length(segrows) > 1
                for j in 2:length(segrows)
                    insert!(new_edge_rows, e + j - 1, segrows[j])
                end
            end

            new_face_rows = _replace_edge_indices(frows, e, length(segrows))

            Delta0 = _build_Delta0(new_edge_rows, size(W, 2))
            Delta1 = _build_Delta1(new_face_rows, length(new_edge_rows))

            changed = true
            progress = true
            break
        end
    end

    nv = size(W, 2)
    ne = size(Delta0, 1)

    gv = hasproperty(block, :global_vertex_ids) ? collect(block.global_vertex_ids) : collect(1:nv)
    ge = hasproperty(block, :global_edge_ids) ? collect(block.global_edge_ids) : collect(1:ne)

    # if new vertices were inserted, extend global ids with fresh local placeholders
    if length(gv) < nv
        append!(gv, collect((maximum(gv) + 1):(maximum(gv) + (nv - length(gv)))))
    end
    if length(ge) < ne
        append!(ge, collect((maximum(ge) + 1):(maximum(ge) + (ne - length(ge)))))
    end

    return (
        id = hasproperty(block, :id) ? block.id : 0,
        W = W,
        Delta0 = Delta0,
        Delta1 = Delta1,
        global_vertex_ids = gv,
        global_edge_ids = ge,
        changed = changed
    )end

"""
    repair_all_original_edges(blocks, FE_orig, lar; atol)

Repair all pairwise T-junction mismatches along shared original edges.
Returns repaired blocks and a summary report.
"""
function repair_all_original_edges(blocks::Vector,
                                   FE_orig::SparseMatrixCSC,
                                   lar;
                                   atol::Float64=1e-7)

    current = Any[]
    for (i, b) in enumerate(blocks)
        nv = size(b.W, 2)
        ne = size(b.Delta0, 1)

        gv = hasproperty(b, :global_vertex_ids) ? collect(b.global_vertex_ids) : collect(1:nv)
        ge = hasproperty(b, :global_edge_ids) ? collect(b.global_edge_ids) : collect(1:ne)

        push!(current, (
            id = hasproperty(b, :id) ? b.id : i,
            W = Matrix{Float64}(b.W),
            Delta0 = SparseMatrixCSC{Int8,Int64}(b.Delta0),
            Delta1 = SparseMatrixCSC{Int8,Int64}(b.Delta1),
            global_vertex_ids = gv,
            global_edge_ids = ge
        ))
    end
    
    repairs = Any[]
    repaired_edges = Int[]

    for e in 1:size(FE_orig, 2)
        incident = edge_incident_blocks(FE_orig, e)
        if length(incident) != 2
            push!(repairs, (
                edge = e,
                ok = false,
                repaired = false,
                reason = "malformed original edge incidence",
                incident_blocks = incident
            ))
            continue
        end

        σ1, σ2 = incident
        seg = edge_segment_from_lar(lar, e)

        cmp_before = compare_edge_subdivision(current[σ1], current[σ2], seg.p, seg.q; atol=atol)

        if cmp_before.ok
            push!(repairs, (
                edge = e,
                ok = true,
                repaired = false,
                incident_blocks = (σ1, σ2),
                before = cmp_before,
                after = cmp_before
            ))
            continue
        end

        target = _unique_sorted(vcat(cmp_before.tvals1, cmp_before.tvals2); atol=atol)

        new1 = _repair_block_edge_to_tvals(current[σ1], seg.p, seg.q, target; atol=atol)
        new2 = _repair_block_edge_to_tvals(current[σ2], seg.p, seg.q, target; atol=atol)

        current[σ1] = new1
        current[σ2] = new2

        cmp_after = compare_edge_subdivision(current[σ1], current[σ2], seg.p, seg.q; atol=atol)

        if cmp_after.ok
            push!(repaired_edges, e)
        end

        push!(repairs, (
            edge = e,
            ok = cmp_after.ok,
            repaired = true,
            incident_blocks = (σ1, σ2),
            target_tvals = target,
            before = cmp_before,
            after = cmp_after
        ))
    end

    return (
        blocks = current,
        ok = all(r.ok for r in repairs),
        repaired_edges = repaired_edges,
        repairs = repairs
    )
end

"""
    assemble_repaired_blocks(blocks)

Reassemble repaired blocks into a single block-diagonal dataset (W, Delta_0, Delta_1).
"""
function assemble_repaired_blocks(blocks::Vector)
    d = size(blocks[1].W, 1)

    total_nv = sum(size(b.W, 2) for b in blocks)
    total_ne = sum(size(b.Delta0, 1) for b in blocks)
    total_nf = sum(size(b.Delta1, 1) for b in blocks)

    Wbig = Matrix{Float64}(undef, d, total_nv)

    I0 = Int[]
    J0 = Int[]
    V0 = Int8[]

    I1 = Int[]
    J1 = Int[]
    V1 = Int8[]

    voff = 0
    eoff = 0
    foff = 0

    for b in blocks
        nv = size(b.W, 2)
        ne = size(b.Delta0, 1)
        nf = size(b.Delta1, 1)

        Wbig[:, voff+1:voff+nv] = b.W

        r0, c0, v0 = findnz(b.Delta0)
        append!(I0, r0 .+ eoff)
        append!(J0, c0 .+ voff)
        append!(V0, Int8.(v0))

        r1, c1, v1 = findnz(b.Delta1)
        append!(I1, r1 .+ foff)
        append!(J1, c1 .+ eoff)
        append!(V1, Int8.(v1))

        voff += nv
        eoff += ne
        foff += nf
    end

    Delta0big = sparse(I0, J0, V0, total_ne, total_nv)
    Delta1big = sparse(I1, J1, V1, total_nf, total_ne)

    return (W = Wbig, Delta_0 = Delta0big, Delta_1 = Delta1big)
end

"""
    repair_report(summary)

Compact text report for `repair_all_original_edges`.
"""
function repair_report(summary)
    lines = String[]
    push!(lines, "Original-edge repair summary")
    push!(lines, "  all repaired checks passed : $(summary.ok)")
    push!(lines, "  repaired original edges    : $(summary.repaired_edges)")

    for r in summary.repairs
        if get(r, :repaired, false) && !get(r, :ok, false)
            push!(lines, "")
            push!(lines, "  edge $(r.edge): still failing after repair")
            push!(lines, "    incident blocks = $(r.incident_blocks)")
            if haskey(r, :target_tvals)
                push!(lines, "    target t-values = $(r.target_tvals)")
            end
            push!(lines, "    before = $(r.before.reason)")
            push!(lines, "    after  = $(r.after.reason)")
            push!(lines, "    after tvals block1 = $(r.after.tvals1)")
            push!(lines, "    after tvals block2 = $(r.after.tvals2)")
        end
    end

    return join(lines, "\n")
end

end # module
module BlockBoundaries

using SparseArrays

export extract_blocks,
       block_boundary_edges,
       block_boundary_vertices,
       block_boundary_loops,
       validate_block_boundary,
       summarize_block_boundaries,
       boundary_report

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

"""
    _normalize_points(W)

Return a point matrix with columns = points.
Accept either d×n or n×d.
"""
function _normalize_points(W::AbstractMatrix)
    s = size(W)
    return s[1] <= s[2] ? Matrix(W) : Matrix(transpose(W))
end

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

"""
    _connected_components(adj)

Connected components of an undirected graph given by adjacency lists.
"""
function _connected_components(adj::Vector{Vector{Int}})
    n = length(adj)
    visited = falses(n)
    comps = Vector{Vector{Int}}()

    for s in 1:n
        if visited[s]
            continue
        end
        stack = [s]
        visited[s] = true
        comp = Int[]
        while !isempty(stack)
            u = pop!(stack)
            push!(comp, u)
            for v in adj[u]
                if !visited[v]
                    visited[v] = true
                    push!(stack, v)
                end
            end
        end
        sort!(comp)
        push!(comps, comp)
    end
    return comps
end

"""
    _loop_order_from_edges(EV, edge_ids)

Try to reconstruct one loop ordering from a set of local edge ids.
Returns a named tuple.
"""
function _loop_order_from_edges(EV::Vector{Vector{Int}}, edge_ids::Vector{Int})
    if isempty(edge_ids)
        return (ok=false, reason="empty component", ordered_edges=Int[], ordered_vertices=Int[])
    end

    verts = sort(unique(v for e in edge_ids for v in EV[e]))
    vpos = Dict(v => k for (k, v) in enumerate(verts))
    deg = zeros(Int, length(verts))
    vedges = [Int[] for _ in verts]

    for e in edge_ids
        u, v = EV[e]
        pu, pv = vpos[u], vpos[v]
        deg[pu] += 1
        deg[pv] += 1
        push!(vedges[pu], e)
        push!(vedges[pv], e)
    end

    if !all(==(2), deg)
        bad = verts[findall(!=(2), deg)]
        return (ok=false,
                reason="boundary component has vertices of degree != 2",
                ordered_edges=Int[],
                ordered_vertices=Int[],
                bad_vertices=bad)
    end

    start_edge = edge_ids[1]
    u0, v0 = EV[start_edge]

    ordered_edges = Int[start_edge]
    ordered_vertices = Int[u0, v0]
    used = Set([start_edge])

    current_v = v0

    while length(used) < length(edge_ids)
        candidates = [e for e in vedges[vpos[current_v]] if !(e in used)]
        if length(candidates) != 1
            return (ok=false,
                    reason="cannot continue unique loop traversal",
                    ordered_edges=ordered_edges,
                    ordered_vertices=ordered_vertices)
        end
        e = candidates[1]
        a, b = EV[e]
        next_v = (a == current_v) ? b : a
        push!(ordered_edges, e)
        push!(ordered_vertices, next_v)
        push!(used, e)
        current_v = next_v
    end

    if ordered_vertices[end] != ordered_vertices[1]
        return (ok=false,
                reason="component is not closed after traversal",
                ordered_edges=ordered_edges,
                ordered_vertices=ordered_vertices)
    end

    return (ok=true,
            reason="simple closed loop",
            ordered_edges=ordered_edges,
            ordered_vertices=ordered_vertices)
end

# ------------------------------------------------------------
# Local block extraction
# ------------------------------------------------------------

"""
    extract_blocks(W, Delta_0, Delta_1)

Extract contiguous local blocks from a block-diagonal dataset.

Conventions:
- `W` is d×nV or nV×d
- `Delta_0` is nE × nV   (edge-vertex)
- `Delta_1` is nF × nE   (face-edge)

Returns a vector of named tuples, one per original face block.
"""
function extract_blocks(W::AbstractMatrix,
                        Delta_0::SparseMatrixCSC,
                        Delta_1::SparseMatrixCSC;
                        nverts::Vector{Int},
                        nedges::Vector{Int},
                        nfaces::Vector{Int})

    P = _normalize_points(W)   # d × nV

    sum(nverts) == size(P, 2) ||
        error("sum(nverts) = $(sum(nverts)) does not match number of vertices in W = $(size(P,2))")

    sum(nedges) == size(Delta_0, 1) ||
        error("sum(nedges) = $(sum(nedges)) does not match number of edge rows in Delta_0 = $(size(Delta_0,1))")

    sum(nverts) == size(Delta_0, 2) ||
        error("sum(nverts) = $(sum(nverts)) does not match number of vertex cols in Delta_0 = $(size(Delta_0,2))")

    sum(nfaces) == size(Delta_1, 1) ||
        error("sum(nfaces) = $(sum(nfaces)) does not match number of face rows in Delta_1 = $(size(Delta_1,1))")

    sum(nedges) == size(Delta_1, 2) ||
        error("sum(nedges) = $(sum(nedges)) does not match number of edge cols in Delta_1 = $(size(Delta_1,2))")

    nb = length(nverts)
    length(nedges) == nb || error("nverts and nedges must have same length")
    length(nfaces) == nb || error("nverts and nfaces must have same length")

    blocks = NamedTuple[]

    v0 = 1
    e0 = 1
    f0 = 1

    for bid in 1:nb
        nv = nverts[bid]
        ne = nedges[bid]
        nf = nfaces[bid]

        vrange = v0:(v0 + nv - 1)
        erange = e0:(e0 + ne - 1)
        frange = f0:(f0 + nf - 1)

        Wi = P[:, vrange]
        Delta0i = Delta_0[erange, vrange]
        Delta1i = Delta_1[frange, erange]

        push!(blocks, (
            id = bid,
            W = Matrix(Wi),
            Delta0 = SparseMatrixCSC(Delta0i),
            Delta1 = SparseMatrixCSC(Delta1i),
            v_offset = v0,
            e_offset = e0,
            f_offset = f0,
            global_vertex_ids = collect(vrange),
            global_edge_ids = collect(erange),
            global_face_ids = collect(frange),
        ))

        v0 += nv
        e0 += ne
        f0 += nf
    end

    return blocks
end


# ------------------------------------------------------------
# Boundary extraction per block
# ------------------------------------------------------------

"""
    block_boundary_edges(block)

Return local and global ids of boundary edges of the block.

Boundary edges are local edges used by exactly one local arranged face.
"""
function block_boundary_edges(block)
    Delta1 = block.Delta1
    nF, nE = size(Delta1)

    multiplicity = zeros(Int, nE)
    for e in 1:nE
        rows, _ = findnz(Delta1[:, e])
        multiplicity[e] = length(rows)
    end

    local_boundary_edges = findall(==(1), multiplicity)
    global_boundary_edges = block.global_edge_ids[local_boundary_edges]

    return (
        local_boundary_edges = local_boundary_edges,
        global_boundary_edges = global_boundary_edges,
        multiplicity = multiplicity,
    )
end

"""
    block_boundary_vertices(block)

Return local and global ids of vertices incident to boundary edges.
"""
function block_boundary_vertices(block)
    b = block_boundary_edges(block)
    Delta0 = block.Delta0

    verts = sort(unique(v for e in b.local_boundary_edges for v in _row_support(Delta0, e)))
    gverts = block.global_vertex_ids[verts]

    return (
        local_boundary_vertices = verts,
        global_boundary_vertices = gverts,
    )
end

"""
    block_boundary_loops(block)

Try to reconstruct the boundary loops of a block from its boundary edges.
"""
function block_boundary_loops(block)
    Delta0 = block.Delta0
    b = block_boundary_edges(block)

    EV = Vector{Vector{Int}}(undef, size(Delta0, 1))
    for e in 1:size(Delta0, 1)
        EV[e] = sort(_row_support(Delta0, e))
    end

    bedges = b.local_boundary_edges
    if isempty(bedges)
        return (
            ok = false,
            reason = "no boundary edges detected",
            components = Any[],
            loops = Any[],
        )
    end

    bverts = sort(unique(v for e in bedges for v in EV[e]))
    vpos = Dict(v => k for (k, v) in enumerate(bverts))
    adj = [Int[] for _ in bverts]
    for e in bedges
        u, v = EV[e]
        pu, pv = vpos[u], vpos[v]
        push!(adj[pu], pv)
        push!(adj[pv], pu)
    end
    for a in adj
        unique!(a)
        sort!(a)
    end

    comps_pos = _connected_components(adj)
    comp_vertices = [[bverts[k] for k in comp] for comp in comps_pos]

    loops = Any[]
    for comp in comp_vertices
        eset = [e for e in bedges if all(v -> v in comp, EV[e])]
        push!(loops, _loop_order_from_edges(EV, eset))
    end

    ok = all(getfield(L, :ok) for L in loops)

    return (
        ok = ok,
        reason = ok ? "all boundary components are simple closed loops" :
                      "some boundary component is not a simple closed loop",
        components = comp_vertices,
        loops = loops,
        boundary_edges = bedges,
    )
end

# ------------------------------------------------------------
# Boundary validation
# ------------------------------------------------------------

"""
    validate_block_boundary(block)

Validate the boundary support of one local face block.
"""
function validate_block_boundary(block)
    b_edges = block_boundary_edges(block)
    b_verts = block_boundary_vertices(block)
    loops = block_boundary_loops(block)

    odd_global_edges = Int[]
    bad_boundary_global_vertices = Int[]

    for (le, m) in enumerate(b_edges.multiplicity)
        if m == 0 || (m > 1 && isodd(m))
            push!(odd_global_edges, block.global_edge_ids[le])
        end
    end

    if !isempty(loops.boundary_edges)
        Delta0 = block.Delta0
        EV = [sort(_row_support(Delta0, e)) for e in 1:size(Delta0, 1)]
        bverts = b_verts.local_boundary_vertices
        deg = Dict(v => 0 for v in bverts)
        for e in loops.boundary_edges
            u, v = EV[e]
            deg[u] = get(deg, u, 0) + 1
            deg[v] = get(deg, v, 0) + 1
        end
        badloc = sort([v for (v, d) in deg if d != 2])
        bad_boundary_global_vertices = block.global_vertex_ids[badloc]
    end

    ok = loops.ok && isempty(odd_global_edges) && isempty(bad_boundary_global_vertices)

    return (
        ok = ok,
        block_id = block.id,
        boundary = loops,
        odd_global_edges = odd_global_edges,
        bad_boundary_global_vertices = bad_boundary_global_vertices,
        local_boundary_edges = b_edges.local_boundary_edges,
        global_boundary_edges = b_edges.global_boundary_edges,
        local_boundary_vertices = b_verts.local_boundary_vertices,
        global_boundary_vertices = b_verts.global_boundary_vertices,
        multiplicity = b_edges.multiplicity,
    )
end

# ------------------------------------------------------------
# Batch summaries
# ------------------------------------------------------------

"""
    summarize_block_boundaries(blocks)

Validate all blocks and return a summary.
"""
function summarize_block_boundaries(blocks)
    reports = [validate_block_boundary(b) for b in blocks]
    failed = [r.block_id for r in reports if !r.ok]

    return (
        ok = isempty(failed),
        nblocks = length(blocks),
        failed_blocks = failed,
        reports = reports,
    )
end

"""
    boundary_report(summary)

Compact textual report for all block boundary validations.
"""
function boundary_report(summary)
    lines = String[]
    push!(lines, "Block boundary summary")
    push!(lines, "  number of blocks      : $(summary.nblocks)")
    push!(lines, "  all blocks valid      : $(summary.ok)")
    push!(lines, "  failed blocks         : $(summary.failed_blocks)")

    for rep in summary.reports
        push!(lines, "")
        push!(lines, "  block $(rep.block_id): ok = $(rep.ok)")
        push!(lines, "    local boundary edges    : $(length(rep.local_boundary_edges))")
        push!(lines, "    local boundary vertices : $(length(rep.local_boundary_vertices))")
        push!(lines, "    odd/suspicious edges    : $(rep.odd_global_edges)")
        push!(lines, "    bad boundary vertices   : $(rep.bad_boundary_global_vertices)")
        push!(lines, "    boundary loops ok       : $(rep.boundary.ok)")
        push!(lines, "    boundary reason         : $(rep.boundary.reason)")
    end

    return join(lines, "\n")
end

end # module
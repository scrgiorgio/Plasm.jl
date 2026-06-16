module CongruenceChains

using SparseArrays
using NearestNeighbors
using Plasm: cop2lar

export vcongruence,
       cellcongruence,
       chaincongruence,
       congruence_report

function _normalize_points(W::AbstractMatrix)
    s = size(W)
    if s[1] <= s[2]
        return Matrix(W)
    else
        return Matrix(transpose(W))
    end
end

function _adjacency_components(adj::Vector{Vector{Int}})
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

    sort!(comps; by = c -> c[1])
    return comps
end

function _canonical_cell(cell)
    return Tuple(sort(unique(cell)))
end

function vcongruence(W::AbstractMatrix; epsilon::Real=1e-6)
    P = _normalize_points(W)
    n = size(P, 2)

    kdtree = KDTree(P)

    adj = [Int[] for _ in 1:n]
    for i in 1:n
        near = inrange(kdtree, P[:, i], epsilon)
        for j in near
            if j != i
                push!(adj[i], j)
                push!(adj[j], i)
            end
        end
    end
    for i in 1:n
        unique!(adj[i])
        sort!(adj[i])
    end

    vclasses = _adjacency_components(adj)

    V = hcat([vec(sum(P[:, cls], dims=2) / length(cls)) for cls in vclasses]...)

    vertex_map = zeros(Int, n)
    for (gid, cls) in enumerate(vclasses), lid in cls
        vertex_map[lid] = gid
    end

    return V, vclasses, vertex_map
end

function cellcongruence(Delta::SparseMatrixCSC, lower_classes; dim::Int)
    cellarray = cop2lar(Delta)

    nlower = maximum(vcat(lower_classes...))
    lower_map = Vector{Int}(undef, nlower)
    for (gid, cls) in enumerate(lower_classes), lid in cls
        lower_map[lid] = gid
    end

    classes = Dict{Tuple{Vararg{Int}}, Vector{Int}}()

    for (lcid, cell) in enumerate(cellarray)
        mapped = [lower_map[x] for x in cell]
        key = _canonical_cell(mapped)

        if length(key) > dim
            push!(get!(classes, key, Int[]), lcid)
        end
    end

    ordered_keys = collect(keys(classes))
    sort!(ordered_keys; by = x -> (length(x), collect(x)))

    cells = [collect(k) for k in ordered_keys]
    outclasses = [classes[k] for k in ordered_keys]

    nlocal = length(cellarray)
    cell_map = zeros(Int, nlocal)
    for (gid, cls) in enumerate(outclasses), lid in cls
        cell_map[lid] = gid
    end

    return cells, outclasses, cell_map
end

function chaincongruence(W::AbstractMatrix,
                         Delta_0::SparseMatrixCSC,
                         Delta_1::SparseMatrixCSC;
                         epsilon::Real=1e-6)

    V, vclasses, vertex_map = vcongruence(W; epsilon=epsilon)
    EV, eclasses, edge_map = cellcongruence(Delta_0, vclasses; dim=1)

    # Delta_1 is face-edge (rows=faces, cols=edges), so transpose it
	FE, fclasses, face_map = cellcongruence(Delta_1, eclasses; dim=2)
    return (
        V = V,
        EV = EV,
        FE = FE,
        vclasses = vclasses,
        eclasses = eclasses,
        fclasses = fclasses,
        vertex_map = vertex_map,
        edge_map = edge_map,
        face_map = face_map,
    )
end

function congruence_report(res)
    nv_local = length(res.vertex_map)
    ne_local = length(res.edge_map)
    nf_local = length(res.face_map)

    nv = size(res.V, 2)
    ne = length(res.EV)
    nf = length(res.FE)

    lines = String[]
    println("Congruence report")
    println("  local vertices -> global vertices : $nv_local -> $nv")
    println("  local edges    -> global edges    : $ne_local -> $ne")
    println("  local faces    -> global faces    : $nf_local -> $nf")
    println("  number of vertex classes         : $(length(res.vclasses))")
    println("  number of edge classes           : $(length(res.eclasses))")
    println("  number of face classes           : $(length(res.fclasses))")
    return join(lines, "\n")
end

end # module
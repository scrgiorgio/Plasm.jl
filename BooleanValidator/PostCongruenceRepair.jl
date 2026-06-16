module PostCongruenceRepair

using LinearAlgebra
using SparseArrays

export edge_face_multiplicity,
       bad_edges,
       repair_bad_global_edges,
       repair_report

# ------------------------------------------------------------
# basic helpers
# ------------------------------------------------------------

_point(V::AbstractMatrix, i::Int) = V[:, i]

function _segment_parameter(x::AbstractVector, a::AbstractVector, b::AbstractVector)
    d = b - a
    den = dot(d, d)
    den == 0 && return 0.0
    return dot(x - a, d) / den
end

function _point_on_segment(x::AbstractVector, a::AbstractVector, b::AbstractVector; atol::Float64=1e-7)
    d = b - a
    den = dot(d, d)
    den == 0 && return false, 0.0

    t = dot(x - a, d) / den
    proj = a + t * d
    ok = (-atol <= t <= 1 + atol) && norm(x - proj) <= atol
    return ok, t
end

function _collinear_on_segment(x::AbstractVector, a::AbstractVector, b::AbstractVector; atol::Float64=1e-7)
    ok, t = _point_on_segment(x, a, b; atol=atol)
    return ok, t
end

function _normalize_t(t::Float64; atol::Float64=1e-7)
    abs(t) <= atol && return 0.0
    abs(t - 1.0) <= atol && return 1.0
    return t
end

function _unique_sorted_tvals(vals::Vector{Tuple{Float64,Int}}; atol::Float64=1e-7)
    sort!(vals, by=first)
    out = Tuple{Float64,Int}[]
    for v in vals
        if isempty(out) || abs(v[1] - out[end][1]) > atol
            push!(out, v)
        end
    end
    return out
end

# ------------------------------------------------------------
# incidence and bad edges
# ------------------------------------------------------------

"""
    edge_face_multiplicity(FE, ne)

Return:
- `multiplicity[e]` = number of faces using edge `e`
- `incident_faces[e]` = list of face ids using edge `e`
"""
function edge_face_multiplicity(FE::Vector{Vector{Int}}, ne::Int)
    multiplicity = zeros(Int, ne)
    incident_faces = [Int[] for _ in 1:ne]
    for (f, edges) in enumerate(FE)
        for e in edges
            multiplicity[e] += 1
            push!(incident_faces[e], f)
        end
    end
    return multiplicity, incident_faces
end

"""
    bad_edges(FE, ne)

Edges whose face multiplicity is not exactly 2.
"""
function bad_edges(FE::Vector{Vector{Int}}, ne::Int)
    mult, faces = edge_face_multiplicity(FE, ne)
    bad = [e for e in 1:ne if mult[e] != 2]
    return bad, mult, faces
end

# ------------------------------------------------------------
# path reconstruction along a bad segment
# ------------------------------------------------------------

"""
    _candidate_vertices_on_edge(V, EV, FE, e_bad, incident_faces; atol)

Collect vertices lying on the same geometric support as the bad edge,
restricted to vertices appearing in the incident faces of e_bad.
"""
function _candidate_vertices_on_edge(
    V::AbstractMatrix,
    EV::Vector{Vector{Int}},
    FE::Vector{Vector{Int}},
    e_bad::Int,
    incident_faces::Vector{Int};
    atol::Float64=1e-7
)
    u, v = EV[e_bad]
    a = _point(V, u)
    b = _point(V, v)

    face_edges = unique(vcat([FE[f] for f in incident_faces]...))
    face_vertices = unique(vcat([EV[e] for e in face_edges]...))

    vals = Tuple{Float64,Int}[]
    for w in face_vertices
        x = _point(V, w)
        ok, t = _collinear_on_segment(x, a, b; atol=atol)
        if ok
            push!(vals, (_normalize_t(float(t); atol=atol), w))
        end
    end

    # always force endpoints
    push!(vals, (0.0, u))
    push!(vals, (1.0, v))

    vals = _unique_sorted_tvals(vals; atol=atol)
    return vals
end

"""
    _path_from_candidates(cands)

Given sorted `(t, vid)` candidates, return ordered vertex path.
"""
function _path_from_candidates(cands::Vector{Tuple{Float64,Int}})
    vids = [v for (_, v) in cands]
    # remove accidental repeated consecutive ids
    out = Int[]
    for v in vids
        if isempty(out) || out[end] != v
            push!(out, v)
        end
    end
    return out
end

# ------------------------------------------------------------
# edge replacement in FE / EV
# ------------------------------------------------------------

function _edge_lookup(EV::Vector{Vector{Int}})
    dict = Dict{Tuple{Int,Int},Int}()
    for (e, uv) in enumerate(EV)
        u, v = uv
        dict[(u,v)] = e
        dict[(v,u)] = e
    end
    return dict
end

function _ensure_edge!(EV::Vector{Vector{Int}}, lookup::Dict{Tuple{Int,Int},Int}, u::Int, v::Int)
    if haskey(lookup, (u,v))
        return lookup[(u,v)]
    else
        push!(EV, [u,v])
        e = length(EV)
        lookup[(u,v)] = e
        lookup[(v,u)] = e
        return e
    end
end

"""
    _replace_edge_in_face(edges, e_bad, path_edges)

Replace one occurrence of e_bad in a face boundary by path_edges.
"""
function _replace_edge_in_face(edges::Vector{Int}, e_bad::Int, path_edges::Vector{Int})
    out = Int[]
    replaced = false
    for e in edges
        if e == e_bad && !replaced
            append!(out, path_edges)
            replaced = true
        else
            push!(out, e)
        end
    end
    return out, replaced
end

# ------------------------------------------------------------
# main repair
# ------------------------------------------------------------

"""
    repair_bad_global_edges(V, EV, FE; atol=1e-7)

Repair bad global edges after congruence by replacing each invalid direct edge
with the ordered path of collinear subdividing edges between the same endpoints,
using only vertices from the incident faces.

Returns a named tuple:
- `V, EV, FE` repaired complexes
- `repaired_edges`
- `failed_edges`
- `report`
"""
function repair_bad_global_edges(
    V::AbstractMatrix,
    EV::Vector{Vector{Int}},
    FE::Vector{Vector{Int}};
    atol::Float64=1e-7
)
    V2 = Matrix{Float64}(V)
    EV2 = [collect(e) for e in EV]
    FE2 = [collect(f) for f in FE]

    lookup = _edge_lookup(EV2)

    bad, mult, incident = bad_edges(FE2, length(EV2))

    repaired_edges = Int[]
    failed_edges = Any[]
    details = Any[]

    for e_bad in bad
        faces = incident[e_bad]
        uv = EV2[e_bad]
        u, v = uv

        # only meaningful if edge belongs to at least one face
        if isempty(faces)
            push!(failed_edges, e_bad)
            push!(details, (edge=e_bad, ok=false, reason="edge belongs to no face"))
            continue
        end

        cands = _candidate_vertices_on_edge(V2, EV2, FE2, e_bad, faces; atol=atol)
        path_vertices = _path_from_candidates(cands)

        if length(path_vertices) < 2
            push!(failed_edges, e_bad)
            push!(details, (edge=e_bad, ok=false, reason="no replacement path found"))
            continue
        end

        # if the path is just [u,v], there is no finer replacement
        if length(path_vertices) == 2
            push!(failed_edges, e_bad)
            push!(details, (edge=e_bad, ok=false, reason="no subdividing path found"))
            continue
        end

        path_edges = Int[]
        for i in 1:length(path_vertices)-1
            a = path_vertices[i]
            b = path_vertices[i+1]
            ee = _ensure_edge!(EV2, lookup, a, b)
            push!(path_edges, ee)
        end

        touched = Int[]
        ok = true
        for f in faces
            newf, replaced = _replace_edge_in_face(FE2[f], e_bad, path_edges)
            if !replaced
                ok = false
                break
            end
            FE2[f] = newf
            push!(touched, f)
        end

        if ok
            push!(repaired_edges, e_bad)
            push!(details, (
                edge = e_bad,
                ok = true,
                old_edge = (u, v),
                incident_faces = faces,
                path_vertices = path_vertices,
                path_edges = path_edges
            ))
        else
            push!(failed_edges, e_bad)
            push!(details, (
                edge = e_bad,
                ok = false,
                reason = "failed to replace bad edge in all incident faces",
                incident_faces = faces
            ))
        end
    end

    # recompute bad edges after repair
    bad_after, mult_after, incident_after = bad_edges(FE2, length(EV2))

    report = (
        repaired_edges = repaired_edges,
        failed_edges = failed_edges,
        bad_edges_before = bad,
        bad_edges_after = bad_after,
        multiplicity_after = mult_after,
        incident_faces_after = incident_after,
        details = details
    )

    return (V = V2, EV = EV2, FE = FE2, report = report)
end

# ------------------------------------------------------------
# text report
# ------------------------------------------------------------

function repair_report(rep)
    lines = String[]
    writeln("Post-congruence repair report")
    writeln("  bad edges before repair : $(rep.bad_edges_before)")
    writeln("  repaired edges          : $(rep.repaired_edges)")
    writeln("  failed edges            : $(rep.failed_edges)")
    writeln("  bad edges after repair  : $(rep.bad_edges_after)")

    for d in rep.details
        if get(d, :ok, false)
            writeln("")
            writeln("  edge $(d.edge): repaired")
            writeln("    old edge       = $(d.old_edge)")
            writeln("    incident faces = $(d.incident_faces)")
            writeln("    path vertices  = $(d.path_vertices)")
            writeln("    path edges     = $(d.path_edges)")
        else
            writeln("")
            writeln("  edge $(d.edge): failed")
            writeln("    reason = $(d.reason)")
        end
    end

    return join(lines, "\n")
end

end # module


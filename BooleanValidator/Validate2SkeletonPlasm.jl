module Validate2SkeletonPlasm

using SparseArrays

export edge_face_multiplicity,
       face_edge_incidence,
       validate_2skeleton_plasm,
       validation_report_plasm

# ------------------------------------------------------------
# Build EF from FE
# ------------------------------------------------------------

"""
    face_edge_incidence(FE, ne)

Build sparse EF matrix of size (ne x nf), where:
- rows are global edges
- columns are global faces
- EF[e,f] = ±1 if edge e belongs to face f

Since FE is given only as unsigned edge lists, this routine uses +1 entries.
This is enough for multiplicity/open/nonmanifold diagnostics.
"""
function face_edge_incidence(FE::Vector{Vector{Int}}, ne::Int)
    I = Int[]
    J = Int[]
    V = Int8[]
    for (f, edges) in enumerate(FE)
        for e in edges
            push!(I, e)
            push!(J, f)
            push!(V, Int8(1))
        end
    end
    return sparse(I, J, V, ne, length(FE))
end

"""
    edge_face_multiplicity(FE, ne)

Return:
- multiplicity[e] = number of incident faces of edge e
- incident_faces[e] = list of face ids using e
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

# ------------------------------------------------------------
# Main validator
# ------------------------------------------------------------

"""
    validate_2skeleton_plasm(V, EV, FE)

Validation aligned with Plasm semantics.

Returns a named tuple with three layers:

- `open`:
    edges with multiplicity == 1
    These are the real closure defects.

- `nonmanifold`:
    edges with multiplicity >= 3
    These are informative and may be acceptable in Plasm.

- `strict_surface`:
    old-style strict test for a closed 2-manifold-like surface.
    This is *not* the primary notion of validity for general Plasm 2-skeletons.

Overall:
- `ok_open` means no open edges
- `ok_plasm` means no open edges (recommended criterion)
- `ok_strict_surface` means no open edges and no nonmanifold edges
"""
function validate_2skeleton_plasm(V::AbstractMatrix, EV::Vector{Vector{Int}}, FE::Vector{Vector{Int}})
    ne = length(EV)
    nf = length(FE)

    EF = face_edge_incidence(FE, ne)
    multiplicity, incident_faces = edge_face_multiplicity(FE, ne)

    open_edges = [e for e in 1:ne if multiplicity[e] == 1]
    isolated_edges = [e for e in 1:ne if multiplicity[e] == 0]
    nonmanifold_edges = [e for e in 1:ne if multiplicity[e] >= 3]

    # "strict surface" criterion:
    # exactly two incident faces on every used edge
    used_edges = [e for e in 1:ne if multiplicity[e] > 0]
    strict_bad_edges = [e for e in used_edges if multiplicity[e] != 2]

    open = (
        ok = isempty(open_edges),
        edges = open_edges,
        count = length(open_edges),
        incident_faces = incident_faces,
    )

    nonmanifold = (
        ok = isempty(nonmanifold_edges),
        edges = nonmanifold_edges,
        count = length(nonmanifold_edges),
        incident_faces = incident_faces,
    )

    strict_surface = (
        ok = isempty(strict_bad_edges),
        edges = strict_bad_edges,
        count = length(strict_bad_edges),
        reason = isempty(strict_bad_edges) ?
            "every used edge has exactly two incident faces" :
            "some used edges do not have exactly two incident faces",
        EF = EF,
        multiplicity = multiplicity,
        incident_faces = incident_faces,
    )

    return (
        ok_open = open.ok,
        ok_plasm = open.ok,
        ok_strict_surface = strict_surface.ok,
        open = open,
        nonmanifold = nonmanifold,
        strict_surface = strict_surface,
        isolated_edges = isolated_edges,
        multiplicity = multiplicity,
        incident_faces = incident_faces,
        EF = EF,
        ne = ne,
        nf = nf,
    )
end

# ------------------------------------------------------------
# Pretty report
# ------------------------------------------------------------

function validation_report_plasm(rep)
    lines = String[]
    push!(lines, "2-skeleton validation report (Plasm semantics)")
    push!(lines, "  number of edges             : $(rep.ne)")
    push!(lines, "  number of faces             : $(rep.nf)")
    push!(lines, "  open edges (mult = 1)       : $(rep.open.count)")
    push!(lines, "  nonmanifold edges (mult ≥ 3): $(rep.nonmanifold.count)")
    push!(lines, "  isolated edges (mult = 0)   : $(length(rep.isolated_edges))")
    push!(lines, "  plasm-valid (no open edges) : $(rep.ok_plasm)")
    push!(lines, "  strict surface valid        : $(rep.ok_strict_surface)")

    if !isempty(rep.open.edges)
        push!(lines, "")
        push!(lines, "  open edges:")
        push!(lines, "    $(rep.open.edges)")
    end

    if !isempty(rep.nonmanifold.edges)
        push!(lines, "")
        push!(lines, "  nonmanifold edges:")
        push!(lines, "    $(rep.nonmanifold.edges)")
    end

    if !isempty(rep.isolated_edges)
        push!(lines, "")
        push!(lines, "  isolated edges:")
        push!(lines, "    $(rep.isolated_edges)")
    end

    return join(lines, "\n")
end

end # module

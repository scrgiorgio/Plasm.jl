module BlockComplexPipeline

export BlockComplex, PipelineOptions, PipelineReport,
       collect_and_finalize, validate_and_repair_2skeleton,
       build_face_vertices, run_boolean_ready_pipeline

"""
Container for the final global block complex.

Fields
------
- V  :: Vector          global vertex coordinates
- EV :: Vector          edge-vertex incidence
- FE :: Vector          face-edge incidence
- FV :: Vector          face-vertex incidence
- meta :: Dict          optional diagnostics / provenance
"""
struct BlockComplex{TV,TE,TF,TM}
    V::TV
    EV::TE
    FE::TF
    FV::TF
    meta::TM
end

Base.length(c::BlockComplex) = length(c.V)

"""
Options controlling the final collection / repair stage.

Assumptions
-----------
This module is meant to be placed *after* the independent collection of
embedded 2D arrangements. Each arrangement is assumed to contribute a local
2-skeleton and optionally local block data.
"""
Base.@kwdef struct PipelineOptions
    perform_validation::Bool = true
    perform_repair::Bool = true
    require_plasm_valid::Bool = true
    require_strict_surface::Bool = false
    build_FV_from_FE_EV::Bool = true
    deduplicate_vertices::Bool = false
    vertex_tol::Float64 = 1.0e-9
    verbose::Bool = true
end

"""
Structured report returned by the pipeline.
"""
Base.@kwdef struct PipelineReport
    ok::Bool = false
    stage::Symbol = :init
    message::String = ""
    validation_before = nothing
    validation_after = nothing
    nV_before::Int = 0
    nE_before::Int = 0
    nF_before::Int = 0
    nV_after::Int = 0
    nE_after::Int = 0
    nF_after::Int = 0
    notes::Vector{String} = String[]
end

# -----------------------------------------------------------------------------
# Public entry point
# -----------------------------------------------------------------------------

"""
Run the complete post-arrangement pipeline and return a single global complex.

Input
-----
`W, Delta_0, Delta_1` are the global data accumulated from the independent
collection of embedded 2D arrangements.

Expected semantics
------------------
- `W`       : global vertex coordinates / geometric support
- `Delta_0` : edge-vertex relation (or an equivalent edge-level topology)
- `Delta_1` : face-edge relation (or an equivalent face-level topology)

Output
------
Returns `(complex, report)` where
- `complex :: BlockComplex` contains `V, EV, FE, FV`
- `report  :: PipelineReport` summarizes validation and repair

Notes
-----
This module deliberately separates *adapters* from *core logic*.
If your current code stores `Delta_0` and `Delta_1` in a slightly different
shape, you only need to adjust the adapter methods below.
"""
function run_boolean_ready_pipeline(W, Delta_0, Delta_1; kwargs...)
    opts = PipelineOptions(; kwargs...)
    return collect_and_finalize(W, Delta_0, Delta_1, opts)
end

function collect_and_finalize(W, Delta_0, Delta_1, opts::PipelineOptions=PipelineOptions())
    report = PipelineReport(stage = :collect)

    # 1. Adapt block-chain data to a global 2-skeleton representation.
    V, EV, FE = adapt_block_chain_to_2skeleton(W, Delta_0, Delta_1)

    report = PipelineReport(
        stage = :adapted,
        nV_before = length(V),
        nE_before = length(EV),
        nF_before = length(FE),
        notes = String["Adapted (W, Delta_0, Delta_1) to provisional (V, EV, FE)."]
    )

    # Optional global vertex deduplication if needed by the caller.
    if opts.deduplicate_vertices
        V, EV, FE, dedup_note = deduplicate_global_vertices(V, EV, FE; tol = opts.vertex_tol)
        push!(report.notes, dedup_note)
    end

    # 2. Validate and optionally repair the 2-skeleton.
    validation_before = nothing
    validation_after = nothing

    if opts.perform_validation
        validation_before = validate_2skeleton(V, EV, FE)
        push!(report.notes, "Initial validation completed.")
    end

    if opts.perform_repair
        V, EV, FE, validation_after, repair_notes = validate_and_repair_2skeleton(V, EV, FE; opts = opts)
        append!(report.notes, repair_notes)
    else
        validation_after = validation_before
    end

    # 3. Build FV if requested.
    FV = opts.build_FV_from_FE_EV ? build_face_vertices(EV, FE) : Vector{Vector{Int}}()

    # 4. Final acceptance logic.
    final_ok = true
    if opts.require_plasm_valid && validation_after !== nothing
        final_ok &= getproperty(validation_after, :ok_plasm)
    end
    if opts.require_strict_surface && validation_after !== nothing
        final_ok &= getproperty(validation_after, :ok_strict_surface)
    end

    complex = BlockComplex(V, EV, FE, FV, Dict(
        :validation_before => validation_before,
        :validation_after => validation_after,
        :notes => copy(report.notes),
    ))

    report = PipelineReport(
        ok = final_ok,
        stage = :done,
        message = final_ok ? "Global complex ready for Boolean stage." :
                             "Global complex built, but final acceptance conditions failed.",
        validation_before = validation_before,
        validation_after = validation_after,
        nV_before = report.nV_before,
        nE_before = report.nE_before,
        nF_before = report.nF_before,
        nV_after = length(V),
        nE_after = length(EV),
        nF_after = length(FE),
        notes = report.notes,
    )

    return complex, report
end

# -----------------------------------------------------------------------------
# Validation / repair
# -----------------------------------------------------------------------------

"""
Validate the 2-skeleton using the user's existing validator if available.

The implementation tries `V2P.validate_2skeleton_plasm` first.
Replace or extend this adapter if your validator lives elsewhere.
"""
function validate_2skeleton(V, EV, FE)
    if isdefined(Main, :V2P) && hasproperty(Main.V2P, :validate_2skeleton_plasm)
        return Main.V2P.validate_2skeleton_plasm(V, EV, FE)
    else
        error("No validator found. Expected Main.V2P.validate_2skeleton_plasm(V,EV,FE).")
    end
end

"""
Check and, if necessary, repair a provisional global 2-skeleton.

Repair policy
-------------
- keep Plasm semantics
- remove open edges if a repair routine is available
- preserve strongly nonmanifold edges
- do not force strict 2-manifoldness
"""
function validate_and_repair_2skeleton(V, EV, FE; opts::PipelineOptions=PipelineOptions())
    notes = String[]
    rep0 = validate_2skeleton(V, EV, FE)
    push!(notes, summary_line(rep0, "before repair"))

    if getproperty(rep0, :ok_plasm)
        push!(notes, "No repair needed: provisional global complex is already Plasm-valid.")
        return V, EV, FE, rep0, notes
    end

    # Try a user-provided repair hook first.
    if isdefined(Main, :repair_global_2skeleton)
        V2, EV2, FE2 = Main.repair_global_2skeleton(V, EV, FE)
        rep1 = validate_2skeleton(V2, EV2, FE2)
        push!(notes, "Applied Main.repair_global_2skeleton.")
        push!(notes, summary_line(rep1, "after repair"))
        return V2, EV2, FE2, rep1, notes
    end

    # Fallback: if the existing pipeline exposes a known repair entry point, adapt here.
    # This is intentionally explicit rather than magical.
    error("Complex is not Plasm-valid and no repair hook is defined. Define Main.repair_global_2skeleton(V,EV,FE).")
end

summary_line(rep, label) = string(
    "Validation ", label, ": ok_open=", getproperty(rep, :ok_open),
    ", ok_plasm=", getproperty(rep, :ok_plasm),
    ", ok_strict_surface=", getproperty(rep, :ok_strict_surface)
)

# -----------------------------------------------------------------------------
# Adapters from block-chain data to (V, EV, FE)
# -----------------------------------------------------------------------------

"""
Adapter from `(W, Delta_0, Delta_1)` to `(V, EV, FE)`.

Default assumption
------------------
The global block-chain complex already stores:
- `W`       as vertex coordinates,
- `Delta_0` as edge-vertex incidence,
- `Delta_1` as face-edge incidence.

If your actual internal representation differs, modify only this function.
"""
function adapt_block_chain_to_2skeleton(W, Delta_0, Delta_1)
    V  = W
    EV = Delta_0
    FE = Delta_1
    return V, EV, FE
end

# -----------------------------------------------------------------------------
# FV construction
# -----------------------------------------------------------------------------

"""
Construct face-vertex incidence from edge-vertex and face-edge incidence.

This implementation assumes each face is given by a cyclic chain of edges, but
it does not attempt to geometrically sort an unordered edge set. If your `FE`
contains unordered edge sets, replace `ordered_face_vertices` with your current
cycle reconstruction routine.
"""
function build_face_vertices(EV, FE)
    return [ordered_face_vertices(EV, f_edges) for f_edges in FE]
end

function ordered_face_vertices(EV, f_edges)
    isempty(f_edges) && return Int[]

    # Expand edges as pairs.
    pairs = [Tuple(EV[e]) for e in f_edges]

    # Start with the first edge.
    a, b = pairs[1]
    used = falses(length(pairs))
    used[1] = true
    verts = [a, b]

    # Greedy edge chaining.
    while count(used) < length(pairs)
        tail = verts[end]
        advanced = false
        for i in eachindex(pairs)
            used[i] && continue
            u, v = pairs[i]
            if u == tail
                push!(verts, v)
                used[i] = true
                advanced = true
                break
            elseif v == tail
                push!(verts, u)
                used[i] = true
                advanced = true
                break
            end
        end
        advanced || break
    end

    # Remove duplicated closure vertex if present.
    if length(verts) > 1 && first(verts) == last(verts)
        pop!(verts)
    end
    return verts
end

# -----------------------------------------------------------------------------
# Optional geometric cleanup
# -----------------------------------------------------------------------------

function deduplicate_global_vertices(V, EV, FE; tol=1.0e-9)
    # Conservative default: exact identity unless the user customizes it.
    # This avoids accidental geometric collapse.
    keyfun(p) = p
    map_old_to_new = Dict{Int,Int}()
    uniq = typeof(V)()
    seen = Dict{Any,Int}()

    for i in eachindex(V)
        k = keyfun(V[i])
        if haskey(seen, k)
            map_old_to_new[i] = seen[k]
        else
            push!(uniq, V[i])
            seen[k] = length(uniq)
            map_old_to_new[i] = length(uniq)
        end
    end

    EV2 = [map(v -> map_old_to_new[v], e) for e in EV]
    note = "Global vertex deduplication applied conservatively (exact equality)."
    return uniq, EV2, FE, note
end

# -----------------------------------------------------------------------------
# Pretty reporting helpers
# -----------------------------------------------------------------------------

function Base.show(io::IO, r::PipelineReport)
    println(io, "PipelineReport(")
    println(io, "  ok          = ", r.ok)
    println(io, "  stage       = ", r.stage)
    println(io, "  message     = \"", r.message, "\"")
    println(io, "  sizes before = (V=", r.nV_before, ", E=", r.nE_before, ", F=", r.nF_before, ")")
    println(io, "  sizes after  = (V=", r.nV_after,  ", E=", r.nE_after,  ", F=", r.nF_after,  ")")
    println(io, "  notes       = ", r.notes)
    print(io, ")")
end

end # module

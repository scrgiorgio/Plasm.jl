using DataStructures
using Random

# -------------------------
# Custom ordering by first component
# -------------------------

struct ByFirstComponent <: Base.Order.Ordering end

Base.Order.lt(::ByFirstComponent, a::Tuple, b::Tuple) = a[1] < b[1]  # change > to > for reversed

# -------------------------
# BoundedBinaryHeap definition
# -------------------------

struct BoundedBinaryHeap{T, O <: Base.Order.Ordering} <: DataStructures.AbstractHeap{T}
    ordering::O
    valtree::Vector{T}
    n::Int # maximum length

    function BoundedBinaryHeap{T}(n::Integer, ordering::Base.Ordering) where T
        n ≥ 1 || throw(ArgumentError("max heap size $n must be ≥ 1"))
        new{T, typeof(ordering)}(ordering, sizehint!(Vector{T}(), n), n)
    end

    function BoundedBinaryHeap{T}(n::Integer, ordering::Base.Ordering, xs::AbstractVector) where T
        n ≥ length(xs) || throw(ArgumentError("initial array is larger than max heap size $n"))
        valtree = sizehint!(DataStructures.heapify(xs, ordering), n)
        new{T, typeof(ordering)}(ordering, valtree, n)
    end
end

# Extra constructors for convenience
BoundedBinaryHeap(n::Integer, ordering::Base.Order.Ordering, xs::AbstractVector{T}) where T =
    BoundedBinaryHeap{T}(n, ordering, xs)

BoundedBinaryHeap{T, O}(n::Integer) where {T, O<:Base.Order.Ordering} =
    BoundedBinaryHeap{T}(n, O())

BoundedBinaryHeap{T, O}(n::Integer, xs::AbstractVector) where {T, O<:Base.Order.Ordering} =
    BoundedBinaryHeap{T}(n, O(), xs)

# Add missing constructor (this one caused the error)
BoundedBinaryHeap{T, O}(n::Integer, ordering::O) where {T, O<:Base.Order.Ordering} =
    BoundedBinaryHeap{T}(n, ordering)

# Interface methods
Base.length(h::BoundedBinaryHeap) = length(h.valtree)
Base.isempty(h::BoundedBinaryHeap) = isempty(h.valtree)
@inline Base.first(h::BoundedBinaryHeap) = h.valtree[1]
Base.pop!(h::BoundedBinaryHeap) = DataStructures.heappop!(h.valtree, h.ordering)

function Base.push!(h::BoundedBinaryHeap, v)
    if length(h) < h.n
        DataStructures.heappush!(h.valtree, v, h.ordering)
    elseif Base.Order.lt(h.ordering, @inbounds(h.valtree[1]), v)
        DataStructures.percolate_down!(h.valtree, 1, v, h.ordering)
    end
    return h
end

# -------------------------
# Example usage
# -------------------------

Random.seed!(1234)  # reproducible example

# Create heap with top 5 (value, metadata) pairs
heap = BoundedBinaryHeap{Tuple{Float64, Vector{Float64}}}(8, ByFirstComponent())

for i in 1:100
    val = randn()
    meta = randn(3)
    push!(heap, (val, meta))
end

println("\nTop entries in heap (by largest value):")
for (val, meta) in heap.valtree
    println("value = ", val, ", metadata = ", meta)
end

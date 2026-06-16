module RobustFaceAffine

using LinearAlgebra
using Main.RobustFaceAffine

export poly_to_plane_affine,
       apply_affine,
       face_to_local2d,
       local2d_to_face3d,
       max_plane_residual

"""
    poly_to_plane_affine(V; T=Float64, atol=nothing)

Return `(A, Ainv)` where `A` is a 4x4 affine transform mapping a 3D polygon
to a local frame whose best-fitting plane is `z ≈ 0`, and `Ainv` maps back.

Input:
- `V` is a `3 x n` matrix of polygon vertices (columns are points)

Properties:
- best normal computed by Newell's method
- origin is the first polygon vertex
- local x-axis is chosen robustly in the plane
- local y-axis completes a right-handed orthonormal frame
"""
function poly_to_plane_affine(V::AbstractMatrix; T=Float64, atol=nothing)
    size(V,1) == 3 || throw(ArgumentError("V must be a 3 x n matrix"))
    n = size(V,2)
    n >= 3 || throw(ArgumentError("Need at least 3 vertices"))

    atol === nothing && (atol = sqrt(eps(T)))

    W = Matrix{T}(V)

    # ------------------------------------------------------------
    # 1. Newell normal
    # ------------------------------------------------------------
    nx = zero(T)
    ny = zero(T)
    nz = zero(T)

    @inbounds for i in 1:n
        j = (i == n) ? 1 : i + 1
        px, py, pz = W[1,i], W[2,i], W[3,i]
        qx, qy, qz = W[1,j], W[2,j], W[3,j]

        nx += (py - qy) * (pz + qz)
        ny += (pz - qz) * (px + qx)
        nz += (px - qx) * (py + qy)
    end

    N = T[nx, ny, nz]

    # ------------------------------------------------------------
    # 2. Fallback if Newell is nearly degenerate
    # ------------------------------------------------------------
    if norm(N) <= atol
        p1 = W[:,1]
        best = zero(T)
        bestN = zeros(T, 3)

        for i in 2:n-1
            for j in i+1:n
                a = W[:,i] - p1
                b = W[:,j] - p1
                c = cross(a, b)
                area2 = norm(c)
                if area2 > best
                    best = area2
                    bestN .= c
                end
            end
        end

        best <= atol && throw(ArgumentError("Degenerate polygon: cannot determine a stable normal"))
        N .= bestN
    end

    zhat = N / norm(N)

    # ------------------------------------------------------------
    # 3. Robust in-plane basis
    # ------------------------------------------------------------
    ex = T[1,0,0]
    ey = T[0,1,0]
    ez = T[0,0,1]

    dx = abs(dot(zhat, ex))
    dy = abs(dot(zhat, ey))
    dz = abs(dot(zhat, ez))

    ref = dx <= dy && dx <= dz ? ex : (dy <= dz ? ey : ez)

    xhat = cross(ref, zhat)
    norm(xhat) <= atol && throw(ArgumentError("Failed to build stable in-plane x-axis"))
    xhat ./= norm(xhat)

    yhat = cross(zhat, xhat)
    yhat ./= norm(yhat)

    # ------------------------------------------------------------
    # 4. Affine transform: p_local = R * (p - p0)
    # ------------------------------------------------------------
    R = [xhat';
         yhat';
         zhat']

    p0 = W[:,1]
    t = -R * p0

    A = Matrix{T}(I, 4, 4)
    A[1:3,1:3] .= R
    A[1:3,4]   .= t

    # inverse: p_world = R' * p_local + p0
    Ainv = Matrix{T}(I, 4, 4)
    Ainv[1:3,1:3] .= transpose(R)
    Ainv[1:3,4]   .= p0

    return A, Ainv
end

"""
    apply_affine(A, V)

Apply 4x4 affine transform `A` to a `3 x n` point matrix `V`.
Returns a `3 x n` matrix.
"""
function apply_affine(A::AbstractMatrix, V::AbstractMatrix)
    size(A) == (4,4) || throw(ArgumentError("A must be 4x4"))
    size(V,1) == 3 || throw(ArgumentError("V must be 3 x n"))

    T = promote_type(eltype(A), eltype(V))
    H = vcat(Matrix{T}(V), ones(T, 1, size(V,2)))
    Y = A * H
    return Y[1:3, :]
end

"""
    face_to_local2d(V; T=Float64, atol=nothing)

Map a 3D face `V` (3 x n) to local coordinates.
Returns `(U, A, Ainv)` where:
- `U` is `2 x n` local planar coordinates,
- `A` maps 3D -> local 3D (with z ≈ 0),
- `Ainv` maps local 3D -> world 3D.
"""
function face_to_local2d(V::AbstractMatrix; T::Type=Float64, atol=nothing)
    A, Ainv = poly_to_plane_affine(V; T=T, atol=atol)
    Vloc3 = apply_affine(A, V)
    U = Vloc3[1:2, :]
    return U, A, Ainv
end

"""
    local2d_to_face3d(U, Ainv)

Map local planar coordinates `U` (2 x n) back to 3D world coordinates
using inverse affine transform `Ainv`.
"""
function local2d_to_face3d(U::AbstractMatrix, Ainv::AbstractMatrix)
    size(U,1) == 2 || throw(ArgumentError("U must be 2 x n"))
    Z = zeros(promote_type(eltype(U), eltype(Ainv)), 1, size(U,2))
    Vloc3 = vcat(U, Z)
    return apply_affine(Ainv, Vloc3)
end

"""
    max_plane_residual(V, A)

After mapping `V` with affine transform `A`, return the maximum absolute
z-coordinate. Useful to measure how well the polygon was flattened.
"""
function max_plane_residual(V::AbstractMatrix, A::AbstractMatrix)
    Vloc = apply_affine(A, V)
    return maximum(abs.(Vloc[3, :]))
end

end # module

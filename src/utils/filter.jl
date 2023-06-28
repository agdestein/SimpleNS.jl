"""
    spectral_cutoff(u, M)

Filter `u` (of size `N` times `N`) to `ubar` of size `M` times `M`.
"""
function spectral_cutoff(u, k)
    T = eltype(u)
    N = size(u, 1)
    uhat = fft(u)
    ubarhat = uhat[[1:k; end-k+1:end], [1:k; end-k+1:end]]
    ubar = sqrt(T(2k) / T(N)) .* abs.(ifft(ubarhat))
    ubar
end

"""
    apply_matrix(A, x)

Apply matrix `A` to a collection of vectors `x` of size `(vecsize, d₁, …, dₙ)`.
"""
function apply_matrix(A, x)
    x = Array(x)
    n, s... = size(x)
    x = reshape(x, n, :)
    y = A * x
    reshape(y, size(y, 1), s...)
end

"""
    create_top_hat_pressure(N, M)

`N` fine points and `M` coarse points in each dimension.
"""
function create_top_hat_pressure(N, M)
    s = N ÷ M
    @assert s * M == N

    i = 1:M
    j = reshape(1:M, 1, :)
    ij = @. i + M * (j - 1)
    ij = repeat(ij, 1, 1, s, s)

    k = reshape(1:s, 1, 1, :)
    l = reshape(1:s, 1, 1, 1, :)

    ijkl = @. s * (i - 1) + k + s * N * (j - 1) + N * (l - 1)
    
    z = fill(1 / s^2, N^2)

    sparse(ij[:], ijkl[:], z)
end

"""
    create_top_hat_u(N, M)

`N` fine points and `M` coarse points in each dimension.
"""
function create_top_hat_u(N, M)
    s = N ÷ M
    @assert s * M == N

    i = 1:M
    j = reshape(1:M, 1, :)
    ij = @. i + M * (j - 1)
    ij = repeat(ij, 1, 1, 1, s)

    k = fill(1, 1, 1, 1)
    l = reshape(1:s, 1, 1, 1, :)

    ijkl = @. s * (i - 1) + k + s * N * (j - 1) + N * (l - 1)
    
    z = fill(1 / s, N * M)

    sparse(ij[:], ijkl[:], z, M^2, N^2)
end

"""
    create_top_hat_v(N, M)

`N` fine points and `M` coarse points in each dimension.
"""
function create_top_hat_v(N, M)
    s = N ÷ M
    @assert s * M == N

    i = 1:M
    j = reshape(1:M, 1, :)
    ij = @. i + M * (j - 1)
    ij = repeat(ij, 1, 1, s, 1)

    k = reshape(1:s, 1, 1, :, 1)
    l = fill(1, 1, 1, 1, 1)

    ijkl = @. s * (i - 1) + k + s * N * (j - 1) + N * (l - 1)
    
    z = fill(1 / s, N * M)

    sparse(ij[:], ijkl[:], z, M^2, N^2)
end

"""
    create_top_hat_velocity(N, M)

`N` fine points and `M` coarse points in each dimension.
"""
function create_top_hat_velocity(N, M)
    Wu = create_top_hat_u(N, M)
    Wv = create_top_hat_v(N, M)
    blockdiag(Wu, Wv)
end

"""
    runge_kutta_method(A, b, c, r)

Get Runge Kutta method. The function checks whether the method is explicit.
"""
function runge_kutta_method(A, b, c, r)
    s = size(A, 1)
    s == size(A, 2) == length(b) == length(c) ||
        error("A, b, and c must have the same sizes")
    isexplicit = all(≈(0), UpperTriangular(A))
    isexplicit || error()

    # T = promote_type(eltype(A), eltype(b), eltype(c), typeof(r))
    # TODO: Find where to pass T
    T = Float64
    A = convert(Matrix{T}, A)
    b = convert(Vector{T}, b)
    c = convert(Vector{T}, c)
    r = convert(T, r)

    # Shift Butcher tableau, as A[1, :] is always zero for explicit methods
    A = [A[2:end, :]; b']

    # Vector with time instances (1 is the time level of final step)
    c = [c[2:end]; 1]

    (; A, b, c, r)
end

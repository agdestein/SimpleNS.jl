"""
    convection(model, V, ϕ, setup; bc_vectors, get_jacobian = false)

Evaluate convective terms `c` and, optionally, Jacobian `∇c = ∂c/∂V`, using the convection
model `model`. The convected quantity is `ϕ` (usually `ϕ = V`).

Non-mutating/allocating/out-of-place version.

See also [`convection!`](@ref).
"""
function convection end

function convection(
    ::NoRegConvectionModel,
    V,
    ϕ,
    setup;
    bc_vectors,
    get_jacobian = false,
    newton_factor = false,
)
    # No regularization
    c, ∇c = convection_components(
        V,
        ϕ,
        setup;
        bc_vectors,
        get_jacobian,
        newton_factor,
    )

    c, ∇c
end

"""
    convection!(model, c, ∇c, V, ϕ, setup, cache; bc_vectors, get_jacobian = false)

Evaluate convective terms `c` and, optionally, Jacobian `∇c = ∂c/∂V`, using the convection
model `model`. The convected quantity is `ϕ` (usually `ϕ = V`).

Mutating/non-allocating/in-place version.

See also [`convection`](@ref).
"""
function convection! end

function convection!(
    ::NoRegConvectionModel,
    c,
    ∇c,
    V,
    ϕ,
    setup,
    cache;
    bc_vectors,
    get_jacobian = false,
    newton_factor = false,
)
    (; c3, ∇c3) = cache

    # No regularization
    convection_components!(
        c,
        ∇c,
        V,
        ϕ,
        setup,
        cache;
        bc_vectors,
        get_jacobian,
        newton_factor,
    )

    c, ∇c
end

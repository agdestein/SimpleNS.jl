"""
    convection(model, V, ϕ, setup; bc_vectors, get_jacobian = false)

Evaluate convective terms `c` and, optionally, Jacobian `∇c = ∂c/∂V`, using the convection
model `model`. The convected quantity is `ϕ` (usually `ϕ = V`).
"""
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

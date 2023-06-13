"""
    convection(model, V, ϕ, setup; get_jacobian = false)

Evaluate convective terms `c` and, optionally, Jacobian `∇c = ∂c/∂V`, using the convection
model `model`. The convected quantity is `ϕ` (usually `ϕ = V`).
"""
function convection(
    V,
    ϕ,
    setup;
    get_jacobian = false,
    newton_factor = false,
)
    # No regularization
    c, ∇c = convection_components(
        V,
        ϕ,
        setup;
        get_jacobian,
        newton_factor,
    )

    c, ∇c
end

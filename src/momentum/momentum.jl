"""
    momentum(
        V, ϕ, p, t, setup;
        get_jacobian = false,
        nopressure = false,
        newton_factor = false,
    )

Calculate RHS of momentum equations and, optionally, Jacobian with respect to velocity field.

  - `V`: velocity field
  - `ϕ`: convected field: e.g. ``\\frac{\\partial (\\phi_x V)}{\\partial x} + \\frac{\\partial (\\phi_y V)}{\\partial y}``; usually `ϕ = V` (so `ϕx = u`, `ϕy = v`)
  - `p`: pressure
  - `get_jacobian`: return `∇F = ∂F/∂V`
  - `nopressure`: exclude pressure gradient; in this case input argument `p` is not used
  - `newton_factor`
"""
function momentum(
    V,
    ϕ,
    p,
    t,
    setup;
    get_jacobian = false,
    nopressure = false,
    newton_factor = false,
)
    (; viscosity_model, force, operators) = setup
    (; G) = operators

    # Convection
    c, ∇c = convection(V, ϕ, setup; get_jacobian, newton_factor)

    # Diffusion
    d, ∇d = diffusion(viscosity_model, V, setup; get_jacobian)

    # Body force
    b = bodyforce(force, t, setup)

    # Residual in Finite Volume form, including the pressure contribution
    F = @. -c + d + b

    # Nopressure = false is the most common situation, in which we return the entire
    # right-hand side vector
    if !nopressure
        F = F .- (G * p)
    end

    if get_jacobian
        # Jacobian requested
        # We return only the Jacobian with respect to V (not p)
        ∇F = @. -∇c + ∇d
    else
        ∇F = nothing
    end

    F, ∇F
end

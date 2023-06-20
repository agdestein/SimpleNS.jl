"""
    momentum(
        V, ϕ, p, setup;
        nopressure = false,
    )

Calculate RHS of momentum equations and, optionally, Jacobian with respect to velocity field.

  - `V`: velocity field
  - `ϕ`: convected field: e.g. ``\\frac{\\partial (\\phi_x V)}{\\partial x} + \\frac{\\partial (\\phi_y V)}{\\partial y}``; usually `ϕ = V` (so `ϕx = u`, `ϕy = v`)
  - `p`: pressure
  - `nopressure`: exclude pressure gradient; in this case input argument `p` is not used
"""
function momentum(
    V,
    ϕ,
    p,
    setup;
    nopressure = false,
)
    (; operators, force) = setup
    (; G) = operators

    # Convection
    c = convection(V, ϕ, setup)

    # Diffusion
    d = diffusion(V, setup)

    # Body force
    b = force

    # Residual in Finite Volume form, including the pressure contribution
    F = @. -c + d + b

    # Nopressure = false is the most common situation, in which we return the entire
    # right-hand side vector
    if !nopressure
        F = F .- (G * p)
    end

    F
end

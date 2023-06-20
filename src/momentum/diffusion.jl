"""
    diffusion(V, setup)

Evaluate diffusive terms.
"""
function diffusion(V, setup)
    (; operators, viscosity) = setup
    (; Diff) = operators

    viscosity .* (Diff * V)
end

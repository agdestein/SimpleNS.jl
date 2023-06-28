"""
    compute_conservation(V, setup)

Compute mass, momentum and energy conservation properties of velocity field.
"""
function compute_conservation(V, setup)
    (; grid, operators) = setup
    (; indu, indv, Ω) = grid
    (; M) = operators

    uₕ = V[indu]
    vₕ = V[indv]

    Ωu = Ω[indu]
    Ωv = Ω[indv]

    # Check if new velocity field is divergence free (mass conservation)
    maxdiv = maximum(abs.(M * V))

    # Calculate total momentum
    umom = sum(Ωu .* uₕ)
    vmom = sum(Ωv .* vₕ)

    # Calculate total kinetic energy
    k = 1 / 2 * sum(Ω .* V .^ 2)

    maxdiv, umom, vmom, k
end

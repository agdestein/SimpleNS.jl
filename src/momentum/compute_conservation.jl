"""
    compute_conservation(V, t, setup)

Compute mass, momentum and energy conservation properties of velocity field.
"""
function compute_conservation(V, t, setup)
    (; grid, operators) = setup
    (; indu, indv, Ω, x, y, xp, yp, hx, hy, gx, gy) = grid
    (; M) = operators

    uₕ = @view V[indu]
    vₕ = @view V[indv]

    Ωu = @view Ω[indu]
    Ωv = @view Ω[indv]

    # Check if new velocity field is divergence free (mass conservation)
    maxdiv = maximum(abs.(M * V))

    # Calculate total momentum
    umom = sum(Ωu .* uₕ)
    vmom = sum(Ωv .* vₕ)

    # Calculate total kinetic energy
    k = 1 / 2 * sum(Ω .* V .^ 2)

    maxdiv, umom, vmom, k
end

"""
    operator_regularization(grid)

Build regularization matrices.
"""
function operator_regularization(grid, operators)
    (; Ω, indu, indv) = grid
    (; Dux, Duy, Dvx, Dvy) = operators
    (; Su_ux, Su_uy, Sv_vx, Sv_vy) = operators

    Δ = max_size(grid)
    α = 1 / 16 * Δ^2

    Ωu⁻¹ = 1 ./ Ω[indu]
    Ωv⁻¹ = 1 ./ Ω[indv]

    # Diffusive matrices in finite-difference setting, without viscosity
    Diffu_f = Diagonal(Ωu⁻¹) * (Dux * Su_ux + Duy * Su_uy)
    Diffv_f = Diagonal(Ωv⁻¹) * (Dvx * Sv_vx + Dvy * Sv_vy)

    (; Diffu_f, Diffv_f, α)
end

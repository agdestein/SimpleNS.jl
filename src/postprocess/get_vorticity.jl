"""
    get_vorticity(V, t, setup)

Get vorticity from velocity field.
"""
function get_vorticity(V, t, setup)
    (; grid) = setup
    (; Nx, Ny) = grid

    Nωx = Nx + 1
    Nωy = Ny + 1

    ω = zeros(Nωx, Nωy)

    vorticity!(ω, V, t, setup)
end

"""
    vorticity!(ω, V, t, setup)

Compute vorticity values at pressure midpoints.
This should be consistent with `operator_postprocessing.jl`.
"""
function vorticity!(ω, V, t, setup)
    (; grid, operators) = setup
    (; indu, indv, Nux_in, Nvy_in, Nx, Ny) = grid
    (; Wv_vx, Wu_uy) = operators

    uₕ = @view V[indu]
    vₕ = @view V[indv]
    ω_flat = reshape(ω, length(ω))

    uₕ_in = uₕ
    vₕ_in = vₕ

    # ω_flat .= Wv_vx * vₕ_in - Wu_uy * uₕ_in
    mul!(ω_flat, Wv_vx, vₕ_in) # a = b * c
    mul!(ω_flat, Wu_uy, uₕ_in, -1, 1) # a = -b * c + a

    ω
end

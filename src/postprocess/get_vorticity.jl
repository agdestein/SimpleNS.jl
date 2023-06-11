"""
    get_vorticity(V, t, setup)

Get vorticity from velocity field.
"""
function get_vorticity(V, t, setup)
    (; grid, boundary_conditions) = setup
    (; Nx, Ny) = grid

    if all(==(:periodic), (boundary_conditions.u.x[1], boundary_conditions.v.y[1]))
        Nωx = Nx + 1
        Nωy = Ny + 1
    else
        Nωx = Nx - 1
        Nωy = Ny - 1
    end

    ω = zeros(Nωx, Nωy)

    vorticity!(ω, V, t, setup)
end

"""
    vorticity!(ω, V, t, setup)

Compute vorticity values at pressure midpoints.
This should be consistent with `operator_postprocessing.jl`.
"""
function vorticity!(ω, V, t, setup)
    (; grid, operators, boundary_conditions) = setup
    (; indu, indv, Nux_in, Nvy_in, Nx, Ny) = grid
    (; Wv_vx, Wu_uy) = operators

    uₕ = @view V[indu]
    vₕ = @view V[indv]
    ω_flat = reshape(ω, length(ω))

    if boundary_conditions.u.x[1] == :periodic && boundary_conditions.v.y[1] == :periodic
        uₕ_in = uₕ
        vₕ_in = vₕ
    else
        # Velocity at inner points
        diagpos = 0
        boundary_conditions.u.x[1] == :pressure && (diagpos = 1)
        boundary_conditions.u.x[1] == :periodic && (diagpos = 1)
        B1D = spdiagm(Nx - 1, Nux_in, diagpos => ones(Nx - 1))
        B2D = I(Ny) ⊗ B1D
        uₕ_in = B2D * uₕ

        diagpos = 0
        boundary_conditions.v.y[1] == :pressure && (diagpos = 1)
        boundary_conditions.v.y[1] == :periodic && (diagpos = 1)
        B1D = spdiagm(Ny - 1, Nvy_in, diagpos => ones(Ny - 1))
        B2D = B1D ⊗ I(Nx)
        vₕ_in = B2D * vₕ
    end

    # ω_flat .= Wv_vx * vₕ_in - Wu_uy * uₕ_in
    mul!(ω_flat, Wv_vx, vₕ_in) # a = b * c
    mul!(ω_flat, Wu_uy, uₕ_in, -1, 1) # a = -b * c + a

    ω
end

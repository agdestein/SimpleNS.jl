"""
    operator_interpolation(grid, boundary_conditions)

Construct averaging operators.
"""
function operator_interpolation(grid, boundary_conditions)
    (; Nx, Ny) = grid
    (; Nux_in, Nux_b, Nux_t, Nuy_in, Nuy_b, Nuy_t) = grid
    (; Nvx_in, Nvx_b, Nvx_t, Nvy_in, Nvy_b, Nvy_t) = grid
    (; hx, hy, hxi, hyi) = grid
    (; Buvy, Bvux) = grid

    weight = 1 / 2

    mat_hx = Diagonal(hxi)
    mat_hy = Diagonal(hyi)

    # Periodic boundary conditions
    if boundary_conditions.u.x == (:periodic, :periodic)
        mat_hx2 = spdiagm(Nx + 2, Nx + 2, [hx[end]; hx; hx[1]])
    else
        mat_hx2 = spdiagm(Nx + 2, Nx + 2, [hx[1]; hx; hx[end]])
    end

    if boundary_conditions.v.y == (:periodic, :periodic)
        mat_hy2 = spdiagm(Ny + 2, Ny + 2, [hy[end]; hy; hy[1]])
    else
        mat_hy2 = spdiagm(Ny + 2, Ny + 2, [hy[1]; hy; hy[end]])
    end

    ## Interpolation operators, u-component
    ## Iu_ux
    diag1 = fill(weight, Nux_t - 1)
    I1D = spdiagm(Nux_t - 1, Nux_t, 0 => diag1, 1 => diag1)

    # Boundary conditions
    Iu_ux_bc = bc_general(
        Nux_t,
        Nux_in,
        Nux_b,
        boundary_conditions.u.x[1],
        boundary_conditions.u.x[2],
        hx[1],
        hx[end],
    )

    # Extend to 2D
    Iu_ux = mat_hy ⊗ (I1D * Iu_ux_bc.B1D)
    Iu_ux_bc = (; Iu_ux_bc..., Bbc = mat_hy ⊗ (I1D * Iu_ux_bc.Btemp))

    ## Iv_uy
    diag1 = fill(weight, Nvx_t - 1)
    I1D = spdiagm(Nvx_t - 1, Nvx_t, 0 => diag1, 1 => diag1)
    boundary_conditions.u.x[1] == :pressure && (I1D[1, :] ./= 2)
    boundary_conditions.u.x[2] == :pressure && (I1D[end, :] ./= 2)

    # The restriction is essentially 1D so it can be directly applied to I1D
    I1D = Bvux * I1D * mat_hx2
    I2D = I(Nuy_t - 1) ⊗ I1D

    # Boundary conditions low/up
    Nb = Nuy_in + 1 - Nvy_in
    Iv_uy_bc_lu = bc_general(
        Nuy_in + 1,
        Nvy_in,
        Nb,
        boundary_conditions.v.y[1],
        boundary_conditions.v.y[2],
        hy[1],
        hy[end],
    )
    Iv_uy_bc_lu = (; Iv_uy_bc_lu..., B2D = Iv_uy_bc_lu.B1D ⊗ I(Nvx_in))
    Iv_uy_bc_lu = (; Iv_uy_bc_lu..., Bbc = Iv_uy_bc_lu.Btemp ⊗ I(Nvx_in))

    # Boundary conditions left/right
    Iv_uy_bc_lr = bc_general_stag(
        Nvx_t,
        Nvx_in,
        Nvx_b,
        boundary_conditions.v.x[1],
        boundary_conditions.v.x[2],
        hx[1],
        hx[end],
    )

    # Take I2D into left/right operators for convenience
    Iv_uy_bc_lr = (; Iv_uy_bc_lr..., B2D = I2D * (I(Nuy_t - 1) ⊗ Iv_uy_bc_lr.B1D))
    Iv_uy_bc_lr = (; Iv_uy_bc_lr..., Bbc = I2D * (I(Nuy_t - 1) ⊗ Iv_uy_bc_lr.Btemp))

    # Resulting operator:
    Iv_uy = Iv_uy_bc_lr.B2D * Iv_uy_bc_lu.B2D

    ## Interpolation operators, v-component

    ## Iu_vx
    diag1 = fill(weight, Nuy_t - 1)
    I1D = spdiagm(Nuy_t - 1, Nuy_t, 0 => diag1, 1 => diag1)
    I1D = Buvy * I1D * mat_hy2
    I2D = I1D ⊗ I(Nvx_t - 1)
    boundary_conditions.v.y[1] == :pressure && (I1D[1, :] ./= 2)
    boundary_conditions.v.y[2] == :pressure && (I1D[end, :] ./= 2)

    # Boundary conditions low/up
    Iu_vx_bc_lu = bc_general_stag(
        Nuy_t,
        Nuy_in,
        Nuy_b,
        boundary_conditions.u.y[1],
        boundary_conditions.u.y[2],
        hy[1],
        hy[end],
    )
    Iu_vx_bc_lu = (; Iu_vx_bc_lu..., B2D = I2D * (Iu_vx_bc_lu.B1D ⊗ I(Nvx_t - 1)))
    Iu_vx_bc_lu = (; Iu_vx_bc_lu..., Bbc = I2D * (Iu_vx_bc_lu.Btemp ⊗ I(Nvx_t - 1)))

    # Boundary conditions left/right
    Nb = Nvx_in + 1 - Nux_in
    Iu_vx_bc_lr = bc_general(
        Nvx_in + 1,
        Nux_in,
        Nb,
        boundary_conditions.u.x[1],
        boundary_conditions.u.x[2],
        hx[1],
        hx[end],
    )

    Iu_vx_bc_lr = (; Iu_vx_bc_lr..., B2D = I(Nuy_in) ⊗ Iu_vx_bc_lr.B1D)
    Iu_vx_bc_lr = (; Iu_vx_bc_lr..., Bbc = I(Nuy_in) ⊗ Iu_vx_bc_lr.Btemp)

    # Resulting operator:
    Iu_vx = Iu_vx_bc_lu.B2D * Iu_vx_bc_lr.B2D

    ## Iv_vy
    diag1 = fill(weight, Nvy_t - 1)
    I1D = spdiagm(Nvy_t - 1, Nvy_t, 0 => diag1, 1 => diag1)

    # Boundary conditions
    Iv_vy_bc = bc_general(
        Nvy_t,
        Nvy_in,
        Nvy_b,
        boundary_conditions.v.y[1],
        boundary_conditions.v.y[2],
        hy[1],
        hy[end],
    )

    # Extend to 2D
    Iv_vy = (I1D * Iv_vy_bc.B1D) ⊗ mat_hx
    Iv_vy_bc = (; Iv_vy_bc..., Bbc = (I1D * Iv_vy_bc.Btemp) ⊗ mat_hx)

    ## Group operators
    operators = (;
        Iu_ux,
        Iv_uy,
        Iu_vx,
        Iv_vy,
        Iu_ux_bc,
        Iv_vy_bc,
        Iv_uy_bc_lr,
        Iv_uy_bc_lu,
        Iu_vx_bc_lr,
        Iu_vx_bc_lu,
    )

    operators
end

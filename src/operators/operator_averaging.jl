"""
    operator_averaging(grid, boundary_conditions)

Construct averaging operators.
"""
function operator_averaging(grid, boundary_conditions)
    (; Nux_in, Nux_b, Nux_t, Nuy_in, Nuy_b, Nuy_t) = grid
    (; Nvx_in, Nvx_b, Nvx_t, Nvy_in, Nvy_b, Nvy_t) = grid
    (; hx, hy) = grid

    # Averaging weight:
    weight = 1 / 2

    ## Averaging operators, u-component

    ## Au_ux: evaluate u at ux location
    diag1 = weight * ones(Nux_t - 1)
    A1D = spdiagm(Nux_t - 1, Nux_t, 0 => diag1, 1 => diag1)

    # Boundary conditions
    Au_ux_bc = bc_general(
        Nux_t,
        Nux_in,
        Nux_b,
        boundary_conditions.u.x[1],
        boundary_conditions.u.x[2],
        hx[1],
        hx[end],
    )

    # Extend to 2D
    Au_ux = I(Nuy_in) ⊗ (A1D * Au_ux_bc.B1D)
    Au_ux_bc = (; Au_ux_bc..., Bbc = I(Nuy_in) ⊗ (A1D * Au_ux_bc.Btemp))

    ## Au_uy: evaluate u at uy location
    diag1 = weight * ones(Nuy_t - 1)
    A1D = spdiagm(Nuy_t - 1, Nuy_t, 0 => diag1, 1 => diag1)

    # Boundary conditions
    Au_uy_bc = bc_general_stag(
        Nuy_t,
        Nuy_in,
        Nuy_b,
        boundary_conditions.u.y[1],
        boundary_conditions.u.y[2],
        hy[1],
        hy[end],
    )

    # Extend to 2D
    Au_uy = (A1D * Au_uy_bc.B1D) ⊗ I(Nux_in)
    Au_uy_bc = (; Au_uy_bc..., Bbc = (A1D * Au_uy_bc.Btemp) ⊗ I(Nux_in))

    ## Averaging operators, v-component

    ## Av_vx: evaluate v at vx location
    diag1 = weight * ones(Nvx_t - 1)
    A1D = spdiagm(Nvx_t - 1, Nvx_t, 0 => diag1, 1 => diag1)

    # Boundary conditions
    Av_vx_bc = bc_general_stag(
        Nvx_t,
        Nvx_in,
        Nvx_b,
        boundary_conditions.v.x[1],
        boundary_conditions.v.x[2],
        hx[1],
        hx[end],
    )

    # Extend to 2D
    Av_vx = I(Nvy_in) ⊗ (A1D * Av_vx_bc.B1D)
    Av_vx_bc = (; Av_vx_bc..., Bbc = I(Nvy_in) ⊗ (A1D * Av_vx_bc.Btemp))

    ## Av_vy: evaluate v at vy location
    diag1 = weight * ones(Nvy_t - 1)
    A1D = spdiagm(Nvy_t - 1, Nvy_t, 0 => diag1, 1 => diag1)

    # Boundary conditions
    Av_vy_bc = bc_general(
        Nvy_t,
        Nvy_in,
        Nvy_b,
        boundary_conditions.v.y[1],
        boundary_conditions.v.y[2],
        hy[1],
        hy[end],
    )

    # Extend to 2D
    Av_vy = (A1D * Av_vy_bc.B1D) ⊗ I(Nvx_in)
    Av_vy_bc = (; Av_vy_bc..., Bbc = (A1D * Av_vy_bc.Btemp) ⊗ I(Nvx_in))

    ## Group operators
    operators = (; Au_ux, Au_uy, Av_vx, Av_vy, Au_ux_bc, Au_uy_bc, Av_vx_bc, Av_vy_bc)

    operators
end

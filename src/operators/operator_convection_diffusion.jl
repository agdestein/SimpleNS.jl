"""
    operator_convection_diffusion(grid, viscosity_model)

Construct convection and diffusion operators.
"""
function operator_convection_diffusion(grid, viscosity_model)
    (; Nx, Ny) = grid
    (; Nux_in, Nux_b, Nux_t, Nuy_in, Nuy_b, Nuy_t) = grid
    (; Nvx_in, Nvx_b, Nvx_t, Nvy_in, Nvy_b, Nvy_t) = grid
    (; hxi, hyi, hxd, hyd) = grid
    (; gxi, gyi, gxd, gyd) = grid
    (; Buvy, Bvux) = grid
    (; Re) = viscosity_model

    ## Convection (differencing) operator Cu

    # Calculates difference from pressure points to velocity points
    diag1 = ones(Nux_t - 2)
    D1D = spdiagm(Nux_t - 2, Nux_t - 1, 0 => -diag1, 1 => diag1)
    Cux = I(Nuy_in) ⊗ D1D
    Dux = Diagonal(hyi) ⊗ D1D

    # Calculates difference from corner points to velocity points
    diag1 = ones(Nuy_t - 2)
    D1D = spdiagm(Nuy_t - 2, Nuy_t - 1, 0 => -diag1, 1 => diag1)
    Cuy = D1D ⊗ I(Nux_in)
    Duy = D1D ⊗ Diagonal(gxi)

    # Cu = [Cux Cuy]
    # Du = [Dux Duy]

    ## Convection (differencing) operator Cv

    # Calculates difference from pressure points to velocity points
    diag1 = ones(Nvx_t - 2)
    D1D = spdiagm(Nvx_t - 2, Nvx_t - 1, 0 => -diag1, 1 => diag1)
    Cvx = I(Nvy_in) ⊗ D1D
    Dvx = Diagonal(gyi) ⊗ D1D

    # Calculates difference from corner points to velocity points
    diag1 = ones(Nvy_t - 2)
    D1D = spdiagm(Nvy_t - 2, Nvy_t - 1, 0 => -diag1, 1 => diag1)
    Cvy = D1D ⊗ I(Nvx_in)
    Dvy = D1D ⊗ Diagonal(hxi)

    # Cv = [Cvx Cvy]
    # Dv = [Dvx Dvy]

    ## Diffusion operator (stress tensor), u-component: similar to averaging, but with mesh sizes

    ## Su_ux: evaluate ux
    diag1 = 1 ./ hxd
    S1D = spdiagm(Nux_t - 1, Nux_t, 0 => -diag1, 1 => diag1)

    # Boundary conditions
    Su_ux_bc = bc_general(Nux_t, Nux_in, Nux_b)

    # Extend to 2D
    Su_ux = I(Ny) ⊗ (S1D * Su_ux_bc.B1D)

    ## Su_uy: evaluate uy
    diag1 = 1 ./ gyd
    S1D = spdiagm(Nuy_t - 1, Nuy_t, 0 => -diag1, 1 => diag1)

    # Boundary conditions
    Su_uy_bc = bc_diff_stag(Nuy_t, Nuy_in, Nuy_b)

    # Extend to 2D
    Su_uy = (S1D * Su_uy_bc.B1D) ⊗ I(Nux_in)

    ## Diffusion operator (stress tensor), v-component: similar to averaging!

    ## Sv_vx: evaluate vx
    diag1 = 1 ./ gxd
    S1D = spdiagm(Nvx_t - 1, Nvx_t, 0 => -diag1, 1 => diag1)

    # Boundary conditions
    Sv_vx_bc = bc_diff_stag(Nvx_t, Nvx_in, Nvx_b)

    # Extend to 2D
    Sv_vx = I(Nvy_in) ⊗ (S1D * Sv_vx_bc.B1D)

    ## Sv_vy: evaluate vy
    diag1 = 1 ./ hyd
    S1D = spdiagm(Nvy_t - 1, Nvy_t, 0 => -diag1, 1 => diag1)

    # Boundary conditions
    Sv_vy_bc = bc_general(Nvy_t, Nvy_in, Nvy_b)

    # Extend to 2D
    Sv_vy = (S1D * Sv_vy_bc.B1D) ⊗ I(Nx)

    ## Assemble operators
    Diffu = 1 / Re * (Dux * Su_ux + Duy * Su_uy)
    Diffv = 1 / Re * (Dvx * Sv_vx + Dvy * Sv_vy)
    Diff = blockdiag(Diffu, Diffv)

    ## Group operators
    operators = (; Cux, Cuy, Cvx, Cvy, Diff)

    operators
end

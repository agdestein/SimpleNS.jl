"""
    operator_averaging(grid)

Construct averaging operators.
"""
function operator_averaging(grid)
    T = eltype(grid.x)
    
    (; Nux_in, Nux_b, Nux_t, Nuy_in, Nuy_b, Nuy_t) = grid
    (; Nvx_in, Nvx_b, Nvx_t, Nvy_in, Nvy_b, Nvy_t) = grid

    # Averaging weight:
    weight = T(1 // 2)

    ## Averaging operators, u-component

    ## Au_ux: evaluate u at ux location
    diag1 = fill(weight, Nux_t - 1)
    A1D = spdiagm(Nux_t - 1, Nux_t, 0 => diag1, 1 => diag1)

    # Boundary conditions
    Au_ux_bc = bc_general(T, Nux_t, Nux_in, Nux_b)

    # Extend to 2D
    Au_ux = I(Nuy_in) ⊗ (A1D * Au_ux_bc.B1D)

    ## Au_uy: evaluate u at uy location
    diag1 = fill(weight, Nuy_t - 1)
    A1D = spdiagm(Nuy_t - 1, Nuy_t, 0 => diag1, 1 => diag1)

    # Boundary conditions
    Au_uy_bc = bc_general_stag(T, Nuy_t, Nuy_in, Nuy_b)

    # Extend to 2D
    Au_uy = (A1D * Au_uy_bc.B1D) ⊗ I(Nux_in)

    ## Averaging operators, v-component

    ## Av_vx: evaluate v at vx location
    diag1 = fill(weight, Nvx_t - 1)
    A1D = spdiagm(Nvx_t - 1, Nvx_t, 0 => diag1, 1 => diag1)

    # Boundary conditions
    Av_vx_bc = bc_general_stag(T, Nvx_t, Nvx_in, Nvx_b)

    # Extend to 2D
    Av_vx = I(Nvy_in) ⊗ (A1D * Av_vx_bc.B1D)

    ## Av_vy: evaluate v at vy location
    diag1 = fill(weight, Nvy_t - 1)
    A1D = spdiagm(Nvy_t - 1, Nvy_t, 0 => diag1, 1 => diag1)

    # Boundary conditions
    Av_vy_bc = bc_general(T, Nvy_t, Nvy_in, Nvy_b)

    # Extend to 2D
    Av_vy = (A1D * Av_vy_bc.B1D) ⊗ I(Nvx_in)

    ## Group operators
    operators = (; Au_ux, Au_uy, Av_vx, Av_vy)

    operators
end

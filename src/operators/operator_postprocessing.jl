"""
    operator_postprocessing(grid, boundary_conditions)

Construct postprocessing operators such as vorticity.
"""
function operator_postprocessing(grid, boundary_conditions)
    (; Nx, Ny, gx, gy, gxd, gyd) = grid

    if all(==(:periodic), (boundary_conditions.u.x[1], boundary_conditions.v.y[1]))
        # For entirely periodic BC, covering entire mesh

        # dv/dx, like Sv_vx
        diag = 1 ./ gxd
        ∂x = spdiagm(
            Nx + 1,
            Nx,
            -Nx => diag[[1]],
            -1 => -diag[1:(end-1)],
            0 => diag[1:(end-1)],
            Nx - 1 => -diag[[end - 1]],
        )
        repeat_x = spdiagm(Nx + 1, Nx, -Nx => [1], 0 => ones(Nx))

        # du/dy, like Su_uy
        diag = 1 ./ gyd
        ∂y = spdiagm(
            Ny + 1,
            Ny,
            -Ny => diag[[1]],
            -1 => -diag[1:(end-1)],
            0 => diag[1:(end-1)],
            Ny - 1 => -diag[[end - 1]],
        )
        repeat_y = spdiagm(Ny + 1, Ny, -Ny => [1], 0 => ones(Ny))
    else
        # dv/dx, like Sv_vx
        diag = 1 ./ gy[2:(end-1)]
        ∂y = spdiagm(Ny - 1, Ny, 0 => -diag, 1 => diag)
        repeat_y = I(Ny - 1)

        # du/dy, like Su_uy
        diag = 1 ./ gx[2:(end-1)]
        ∂x = spdiagm(Nx - 1, Nx, 0 => -diag, 1 => diag)
        repeat_x = I(Nx - 1)
    end

    # Extend to 2D
    Wu_uy = ∂y ⊗ repeat_x
    Wv_vx = repeat_y ⊗ ∂x

    (; Wv_vx, Wu_uy)
end

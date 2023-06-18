"""
    operator_postprocessing(grid)

Construct postprocessing operators such as vorticity.
"""
function operator_postprocessing(grid)
    (; Nx, Ny, gxd, gyd) = grid

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

    # Extend to 2D
    Wu_uy = ∂y ⊗ repeat_x
    Wv_vx = repeat_y ⊗ ∂x

    (; Wv_vx, Wu_uy)
end

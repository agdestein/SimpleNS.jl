"""
    get_setup(
        Nx, Ny, Lx, Ly;
        bodyforce_u = (x, y) -> 0.0,
        bodyforce_v = (x, y) -> 0.0,
        viscosity = 1 / 1000,
    )

Create 2D setup.
"""
function get_setup(
    Nx, Ny, xlims, ylims;
    bodyforce_u = (x, y) -> 0,
    bodyforce_v = (x, y) -> 0,
    viscosity = 1 / 1000,
)
    grid = create_grid(Nx, Ny, xlims, ylims)
    operators = get_operators(grid)
    force = create_body_force(bodyforce_u, bodyforce_v, grid)
    (; grid, operators, force, viscosity)
end

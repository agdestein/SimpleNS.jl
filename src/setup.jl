"""
    Setup(grid, viscosity_model, force, operators)
`
Simulation setup.
"""
struct Setup{
    T,
    V<:AbstractViscosityModel{T},
    F<:AbstractBodyForce{T},
}
    grid::Grid{T}
    viscosity_model::V
    force::F
    operators::Operators{T}
end

"""
    Setup(
        Nx, Ny, Lx, Ly;
        viscosity_model = LaminarModel(; Re = 1000.0),
        bodyforce_u = (x, y) -> 0.0,
        bodyforce_v = (x, y) -> 0.0,
        steady_force = true,
    )

Create 2D setup.
"""
function Setup(
    Nx, Ny, xlims, ylims;
    viscosity_model = LaminarModel(; Re = 1000.0),
    bodyforce_u = (x, y) -> 0.0,
    bodyforce_v = (x, y) -> 0.0,
    steady_force = true,
)
    grid = Grid(Nx, Ny, xlims, ylims)
    if steady_force
        force = SteadyBodyForce(bodyforce_u, bodyforce_v, grid)
    else
        force = UnsteadyBodyForce(bodyforce_u, bodyforce_v, grid)
    end
    operators = Operators(grid, viscosity_model)
    Setup(grid, viscosity_model, force, operators)
end

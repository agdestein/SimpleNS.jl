"""
    Setup(grid, boundary_conditions, viscosity_model, force, operators)

Simulation setup.
"""
struct Setup{
    T,
    N,
    B<:BoundaryConditions{T},
    V<:AbstractViscosityModel{T},
    F<:AbstractBodyForce{T},
}
    grid::Grid{T,N}
    boundary_conditions::B
    viscosity_model::V
    force::F
    operators::Operators{T}
end

"""
    Setup(
        x, y;
        viscosity_model = LaminarModel(; Re = 1000.0),
        u_bc = (x, y, t) -> 0.0,
        v_bc = (x, y, t) -> 0.0,
        dudt_bc = nothing,
        dvdt_bc = nothing,
        bc_type = (;
            u = (; x = (:periodic, :periodic), y = (:periodic, :periodic)),
            v = (; x = (:periodic, :periodic), y = (:periodic, :periodic)),
        ),
        bodyforce_u = (x, y) -> 0.0,
        bodyforce_v = (x, y) -> 0.0,
        steady_force = true,
    )

Create 2D setup.
"""
function Setup(
    x,
    y;
    viscosity_model = LaminarModel(; Re = 1000.0),
    u_bc = (x, y, t) -> 0.0,
    v_bc = (x, y, t) -> 0.0,
    dudt_bc = nothing,
    dvdt_bc = nothing,
    bc_type = (;
        u = (; x = (:periodic, :periodic), y = (:periodic, :periodic)),
        v = (; x = (:periodic, :periodic), y = (:periodic, :periodic)),
    ),
    bodyforce_u = (x, y) -> 0.0,
    bodyforce_v = (x, y) -> 0.0,
    steady_force = true,
)
    boundary_conditions =
        BoundaryConditions(u_bc, v_bc; dudt_bc, dvdt_bc, bc_type, T = eltype(x))
    grid = Grid(x, y)
    if steady_force
        force = SteadyBodyForce(bodyforce_u, bodyforce_v, grid)
    else
        force = UnsteadyBodyForce(bodyforce_u, bodyforce_v, grid)
    end
    operators = Operators(grid, boundary_conditions, viscosity_model)
    Setup(grid, boundary_conditions, viscosity_model, force, operators)
end

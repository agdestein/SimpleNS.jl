"""
    BoundaryConditions{T}

Boundary conditions with floating point type `T`.
"""
Base.@kwdef struct BoundaryConditions{T,U,V,W,DU,DV,DW}
    u::NamedTuple = (;)
    v::NamedTuple = (;)
    u_bc::U
    v_bc::V
    dudt_bc::DU
    dvdt_bc::DV
    p_bc::NamedTuple = (;)
end

"""
    BoundaryConditions(u_bc, v_bc; T = Float64, bc_type, kwargs...)

Create discrete boundary condtions.

Values should either be scalars or vectors. All values `(u, v, p, k, e)` are
defined at (x, y) locations, i.e. the corners of pressure volumes, so they
cover the entire domain, including corners.
"""
function BoundaryConditions(
    u_bc,
    v_bc;
    T = Float64,
    bc_type,
    dudt_bc = nothing,
    dvdt_bc = nothing,
    kwargs...,
)
    bc_type.u.x[1] ∈ (:dirichlet, :periodic, :pressure) || error("Wrong BC for u-left")
    bc_type.u.x[2] ∈ (:dirichlet, :periodic, :pressure) || error("Wrong BC for u-right")
    bc_type.u.y[1] ∈ (:dirichlet, :periodic, :symmetric) || error("Wrong BC for u-low")
    bc_type.u.y[2] ∈ (:dirichlet, :periodic, :symmetric) || error("Wrong BC for u-up")

    bc_type.v.x[1] ∈ (:dirichlet, :periodic, :symmetric) || error("Wrong BC for v-left")
    bc_type.v.x[2] ∈ (:dirichlet, :periodic, :symmetric) || error("Wrong BC for v-right")
    bc_type.v.y[1] ∈ (:dirichlet, :periodic, :pressure) || error("Wrong BC for v-low")
    bc_type.v.y[2] ∈ (:dirichlet, :periodic, :pressure) || error("Wrong BC for v-up")

    # Pressure (for boundaries marked with `:pressure`)
    p∞ = zero(T)
    p_bc = (; x = (p∞, p∞), y = (p∞, p∞))

    BoundaryConditions{
        T,
        typeof(u_bc),
        typeof(v_bc),
        Nothing,
        typeof(dudt_bc),
        typeof(dvdt_bc),
        Nothing,
    }(;
        bc_type...,
        u_bc,
        v_bc,
        w_bc = nothing,
        dudt_bc,
        dvdt_bc,
        dwdt_bc = nothing,
        p_bc,
        kwargs...,
    )
end

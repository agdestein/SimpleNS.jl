"""
    create_initial_conditions(
        setup;
        initial_velocity_u,
        initial_velocity_v,
        initial_pressure = nothing,
    )

Create initial vectors. If `p_initial` is a function instead of
`nothing`, calculate compatible IC for the pressure.
"""
function create_initial_conditions(
    setup;
    initial_velocity_u,
    initial_velocity_v,
    initial_pressure = nothing,
)
    (; grid, operators) = setup
    (; xu, yu, xv, yv, xpp, ypp, Ω) = grid
    (; G, M) = operators

    # Allocate velocity and pressure
    u = zero(xu)[:]
    v = zero(xv)[:]
    p = zero(xpp)[:]

    # Initial velocities
    u .= initial_velocity_u.(xu, yu)[:]
    v .= initial_velocity_v.(xv, yv)[:]
    V = [u[:]; v[:]]

    # Kinetic energy and momentum of initial velocity field
    # Iteration 1 corresponds to t₀ = 0 (for unsteady simulations)
    maxdiv, umom, vmom, k = compute_conservation(V, setup)

    if maxdiv > 1e-12
        @warn "Initial velocity field not (discretely) divergence free: $maxdiv.\n" *
              "Performing additional projection."

        # Make velocity field divergence free
        f = M * V
        Δp = pressure_poisson(setup, f)
        V .-= (G * Δp) ./ Ω
    end

    # Initial pressure: should in principle NOT be prescribed (will be calculated if p_initial)
    if isnothing(initial_pressure)
        p = pressure_additional_solve(setup, V, p)
    else
        p .= initial_pressure.(xpp, ypp)[:]
    end

    V, p
end

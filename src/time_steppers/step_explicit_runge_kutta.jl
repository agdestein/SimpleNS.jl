"""
    step(stepper::ExplicitRungeKuttaStepper, Δt; bc_vectors = nothing)

Perform one time step for the general explicit Runge-Kutta method (ERK).

Dirichlet boundary points are not part of solution vector but are prescribed in a strong
manner via the `u_bc` and `v_bc` functions.
"""
function step(stepper::ExplicitRungeKuttaStepper, Δt; bc_vectors = nothing)
    (; method, setup, pressure_solver, n, V, p, t, Vₙ, pₙ, tₙ) = stepper
    (; grid, operators, boundary_conditions) = setup
    (; bc_unsteady) = boundary_conditions
    (; Ω⁻¹) = grid
    (; G, M) = operators
    (; A, b, c, p_add_solve) = method

    # Update current solution (does not depend on previous step size)
    n += 1
    Vₙ = V
    pₙ = p
    tₙ = t
    Δtₙ = Δt

    # Number of stages
    nV = length(V)
    np = length(p)
    nstage = length(b)

    # Reset RK arrays
    tᵢ = tₙ
    kV = zeros(nV, 0)
    kp = zeros(np, 0)

    ## Start looping over stages

    # At i = 1 we calculate F₁, p₂ and u₂
    # ⋮
    # At i = s we calculate Fₛ, pₙ₊₁, and uₙ₊₁
    for i = 1:nstage
        # Right-hand side for tᵢ based on current velocity field uₕ, vₕ at level i. This
        # includes force evaluation at tᵢ and pressure gradient. Boundary conditions will be
        # set through `get_bc_vectors` inside momentum. The pressure p is not important here,
        # it will be removed again in the next step
        F, ∇F = momentum(V, V, p, tᵢ, setup; bc_vectors)

        # Store right-hand side of stage i
        # Remove the -G*p contribution (but not y_p)
        kVᵢ = Ω⁻¹ .* (F + G * p)
        kV = [kV kVᵢ]

        # Update velocity current stage by sum of Fᵢ's until this stage, weighted
        # with Butcher tableau coefficients. This gives uᵢ₊₁, and for i=s gives uᵢ₊₁
        V = kV * A[i, 1:i]

        # Boundary conditions at tᵢ₊₁
        tᵢ = tₙ + c[i] * Δtₙ
        if isnothing(bc_vectors) || bc_unsteady
            bc_vectors = get_bc_vectors(setup, tᵢ)
        end
        (; yM) = bc_vectors

        # Divergence of intermediate velocity field
        f = (M * (Vₙ / Δtₙ + V) + yM / Δtₙ) / c[i]

        # Solve the Poisson equation, but not for the first step if the boundary conditions are steady
        if boundary_conditions.bc_unsteady || i > 1
            p = pressure_poisson(pressure_solver, f)
        else
            # Bc steady AND i = 1
            p = pₙ
        end

        Gp = G * p

        # Update velocity current stage, which is now divergence free
        V = @. Vₙ + Δtₙ * (V - c[i] * Ω⁻¹ * Gp)
    end

    # For steady bc we do an additional pressure solve
    # That saves a pressure solve for i = 1 in the next time step
    if !bc_unsteady || p_add_solve
        p = pressure_additional_solve(pressure_solver, V, p, tₙ + Δtₙ, setup; bc_vectors)
    end

    t = tₙ + Δtₙ

    TimeStepper(; method, setup, pressure_solver, n, V, p, t, Vₙ, pₙ, tₙ)
end

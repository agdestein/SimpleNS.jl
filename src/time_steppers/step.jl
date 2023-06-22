"""
    step(stepper, Δt)

Perform one time step for the general explicit Runge-Kutta method (ERK).
"""
function step(stepper, Δt)
    (; method, setup, n, V, p, t) = stepper
    (; grid, operators) = setup
    (; Ω) = grid
    (; G, M) = operators
    (; A, b, c) = method

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
    kV = repeat(V, 1, 0)
    kp = repeat(p, 1, 0)

    ## Start looping over stages

    # At i = 1 we calculate F₁, p₂ and u₂
    # ⋮
    # At i = s we calculate Fₛ, pₙ₊₁, and uₙ₊₁
    for i = 1:nstage
        # Right-hand side for tᵢ based on current velocity field uₕ, vₕ at level i. This
        # includes force evaluation at tᵢ and pressure gradient. The pressure p is not important here,
        # it will be removed again in the next step
        F = momentum(V, V, p, setup; nopressure = true)

        # Store right-hand side of stage i
        kVᵢ = F ./ Ω
        kV = [kV kVᵢ]

        # Update velocity current stage by sum of Fᵢ's until this stage, weighted
        # with Butcher tableau coefficients. This gives uᵢ₊₁, and for i=s gives uᵢ₊₁
        ΔV = kV * A[i, 1:i]

        # Divergence of intermediate velocity field
        ftemp = @. (Vₙ / Δtₙ + ΔV) / c[i]
        f = M * ftemp

        # Solve the Poisson equation, but not for the first step if the boundary conditions are steady
        if i > 1
            p = pressure_poisson(setup, f)
        else
            # Bc steady AND i = 1
            p = pₙ
        end

        Gp = G * p

        # Update velocity current stage, which is now divergence free
        V = @. Vₙ + Δtₙ * (ΔV - c[i] * Gp / Ω )
    end

    # For steady bc we do an additional pressure solve
    # That saves a pressure solve for i = 1 in the next time step
    p = pressure_additional_solve(setup, V, p)

    t = tₙ + Δtₙ

    (; method, setup, n, V, p, t)
end


"""
    step_rk4(stepper, Δt)

Perform one time step for the general explicit Runge-Kutta method (ERK).
"""
function step_rk4(stepper, Δt)
    (; method, setup, n, V, p, t) = stepper
    (; grid, operators) = setup
    (; Ω) = grid
    (; G, M) = operators

    # Update current solution (does not depend on previous step size)
    n += 1
    Vₙ = V
    pₙ = p
    tₙ = t
    Δtₙ = Δt


    ## Stage 1

    F = momentum(V, V, p, setup; nopressure = true)
    k₁ = F ./ Ω

    ΔV = k₁ / 2
    
    # Divergence of intermediate velocity field
    f = M * ΔV

    # Solve the Poisson equation, but not for the first step if the boundary conditions are steady
    Δp = pₙ

    # Make velocity increment divergence free
    ΔV = ΔV - G * Δp

    # Update velocity current stage, which is now divergence free
    V = @. Vₙ + Δtₙ * ΔV


    ## Stage 2

    F = momentum(V, V, p, setup; nopressure = true)
    k₂ = F ./ Ω

    ΔV = k₂ / 2
    
    # Divergence of intermediate velocity field
    f = M * ΔV

    # Solve the Poisson equation
    Δp = pressure_poisson(setup, f)

    # Make velocity increment divergence free
    ΔV = ΔV - G * Δp

    # Update velocity current stage, which is now divergence free
    V = @. Vₙ + Δtₙ * ΔV


    ## Stage 3

    F = momentum(V, V, p, setup; nopressure = true)
    k₃ = F ./ Ω

    ΔV = k₃
    
    # Divergence of intermediate velocity field
    f = M * ΔV

    # Solve the Poisson equation
    Δp = pressure_poisson(setup, f)

    # Make velocity increment divergence free
    ΔV = ΔV - G * Δp

    # Update velocity current stage, which is now divergence free
    V = @. Vₙ + Δtₙ * ΔV


    ## Stage 4

    F = momentum(V, V, p, setup; nopressure = true)
    k₄ = F ./ Ω

    ΔV = k₁ / 6 + k₂ / 3 + k₃ / 3 + k₄ / 6
    
    # Divergence of intermediate velocity field
    f = M * ΔV

    # Solve the Poisson equation
    Δp = pressure_poisson(setup, f)

    # Make velocity increment divergence free
    ΔV = ΔV - G * Δp

    # Update velocity current stage, which is now divergence free
    V = @. Vₙ + Δtₙ * ΔV


    # For steady bc we do an additional pressure solve
    # That saves a pressure solve for i = 1 in the next time step
    p = pressure_additional_solve(setup, V, p)

    t = tₙ + Δtₙ

    (; method, setup, n, V, p, t)
end

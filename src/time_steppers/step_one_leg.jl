"""
    step(stepper::OneLegStepper, Δt)

Do one time step using one-leg-β-method.
"""
function step(stepper::OneLegStepper, Δt)
    (; method, setup, pressure_solver, n, V, p, t, Vₙ, pₙ, tₙ) = stepper
    (; p_add_solve, β) = method
    (; grid, operators) = setup
    (; G, M) = operators
    (; Ω) = grid

    # Update current solution (does not depend on previous step size)
    Δtₙ₋₁ = t - tₙ
    n += 1
    Vₙ₋₁ = Vₙ
    pₙ₋₁ = pₙ
    Vₙ = V
    pₙ = p
    tₙ = t
    Δtₙ = Δt
    @assert Δtₙ ≈ Δtₙ₋₁

    # Intermediate ("offstep") velocities
    t = tₙ + β * Δtₙ
    V = @. (1 + β) * Vₙ - β * Vₙ₋₁
    p = @. (1 + β) * pₙ - β * pₙ₋₁

    # Right-hand side of the momentum equation
    F, = momentum(V, V, p, t, setup)

    # Take a time step with this right-hand side, this gives an intermediate velocity field
    # (not divergence free)
    V = @. (2β * Vₙ - (β - 1 // 2) * Vₙ₋₁ + Δtₙ  * F / Ω) / (β + 1 // 2)

    # Adapt time step for pressure calculation
    Δtᵦ = Δtₙ / (β + 1 // 2)

    # Divergence of intermediate velocity field
    f = (M * V) / Δtᵦ

    # Solve the Poisson equation for the pressure
    Δp = pressure_poisson(pressure_solver, f)
    GΔp = G * Δp

    # Update velocity field
    V = @. V - Δtᵦ * GΔp / Ω

    # Update pressure (second order)
    p = @. 2pₙ - pₙ₋₁ + 4 // 3 * Δp

    # Alternatively, do an additional Poisson solve
    if p_add_solve
        p = pressure_additional_solve(pressure_solver, V, p, tₙ + Δtₙ, setup)
    end

    t = tₙ + Δtₙ

    TimeStepper(; method, setup, pressure_solver, n, V, p, t, Vₙ, pₙ, tₙ)
end

"""
    step(stepper::AdamsBashforthCrankNicolsonStepper, Δt)

Perform one time step with Adams-Bashforth for convection and Crank-Nicolson for diffusion.

Output includes convection terms at `tₙ`, which will be used in next time step
in the Adams-Bashforth part of the method.
"""
function step(stepper::AdamsBashforthCrankNicolsonStepper, Δt)
    (; method, setup, pressure_solver, n, V, p, t, Vₙ, pₙ, tₙ) = stepper
    (; viscosity_model, force, grid, operators) = setup
    (; NV, Ω) = grid
    (; G, M) = operators
    (; Diff) = operators
    (; p_add_solve, α₁, α₂, θ) = method

    # For the first time step, this might be necessary
    cₙ, = convection(Vₙ, Vₙ, setup)

    # Advance one step
    Δtₙ₋₁ = t - tₙ
    n += 1
    Vₙ = V
    pₙ = p
    tₙ = t
    Δtₙ = Δt
    cₙ₋₁ = cₙ
    @assert Δtₙ ≈ Δtₙ₋₁

    # Evaluate boundary conditions and force at starting point
    bₙ = bodyforce(force, tₙ, setup)

    # Convection of current solution
    cₙ, = convection(Vₙ, Vₙ, setup)

    bₙ₊₁ = bodyforce(force, tₙ + Δt, setup)

    # Crank-Nicolson weighting for force and diffusion boundary conditions
    b = @. (1 - θ) * bₙ + θ * bₙ₊₁

    Gpₙ = G * pₙ

    d = Diff * V

    # Right hand side of the momentum equation update
    Rr = @. Vₙ + 1 / Ω * Δt * (-(α₁ * cₙ + α₂ * cₙ₋₁) + (1 - θ) * d + b - Gpₙ)

    # Implicit time-stepping for diffusion
    if viscosity_model isa LaminarModel
        # Use precomputed LU decomposition
        V = Diff_fact \ Rr
    else
        # Get `∇d` since `Diff` is not constant
        d, ∇d = diffusion(V, t, setup; get_jacobian = true)
        V = ∇d \ Rr
    end

    # Divergence of `Ru` and `Rv` is directly calculated with `M`
    f = (M * V) / Δt

    # Solve the Poisson equation for the pressure
    Δp = pressure_poisson(pressure_solver, f)

    # Update velocity field
    V -= Δt .* 1 ./ Ω .* (G * Δp .+ y_Δp)

    # First order pressure:
    p = pₙ .+ Δp

    if p_add_solve
        p = pressure_additional_solve(pressure_solver, V, p, tₙ + Δt, setup)
    end

    t = tₙ + Δtₙ

    TimeStepper(; method, setup, pressure_solver, n, V, p, t, Vₙ, pₙ, tₙ)
end

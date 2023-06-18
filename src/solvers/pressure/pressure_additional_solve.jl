"""
    pressure_additional_solve(pressure_solver, V, p, t, setup)

Do additional pressure solve. This makes the pressure compatible with the velocity
field, resulting in same order pressure as velocity.
"""
function pressure_additional_solve(pressure_solver, V, p, t, setup)
    (; grid, operators) = setup
    (; Ω) = grid
    (; M) = operators

    # Momentum already contains G*p with the current p, we therefore effectively solve for
    # the pressure difference
    F, = momentum(V, V, p, t, setup)
    f = M * (F ./ Ω)

    Δp = pressure_poisson(pressure_solver, f)
    p .+ Δp
end

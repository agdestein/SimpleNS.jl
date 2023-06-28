"""
    pressure_additional_solve(setup, V, p)

Do additional pressure solve. This makes the pressure compatible with the velocity
field, resulting in same order pressure as velocity.
"""
function pressure_additional_solve(setup, V, p)
    (; grid, operators) = setup;
    (; Ω) = grid;
    (; M) = operators;

    # Momentum already contains G*p with the current p, we therefore effectively solve for
    # the pressure difference
    F = momentum(V, V, p, setup)
    f = M * (F ./ Ω)

    Δp = pressure_poisson(setup, f)
    p .+ Δp
end

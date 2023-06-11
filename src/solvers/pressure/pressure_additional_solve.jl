"""
    pressure_additional_solve(pressure_solver, V, p, t, setup; bc_vectors = nothing)

Do additional pressure solve. This makes the pressure compatible with the velocity
field, resulting in same order pressure as velocity.
"""
function pressure_additional_solve(pressure_solver, V, p, t, setup; bc_vectors = nothing)
    (; grid, operators, boundary_conditions) = setup
    (; bc_unsteady) = boundary_conditions
    (; Ω⁻¹) = grid
    (; M) = operators

    # Get updated BC for ydM (time derivative of BC in ydM)
    # FIXME: `get_bc_vectors` are called to often (also inside `momentum!`)
    if isnothing(bc_vectors) || bc_unsteady
        bc_vectors = get_bc_vectors(setup, t)
    end
    (; ydM) = bc_vectors

    # Momentum already contains G*p with the current p, we therefore effectively solve for
    # the pressure difference
    F, = momentum(V, V, p, t, setup; bc_vectors)
    f = M * (Ω⁻¹ .* F) + ydM

    Δp = pressure_poisson(pressure_solver, f)
    p .+ Δp
end

"""
    get_velocity(V, t, setup)

Get velocity values at pressure points. Interpolate velocities to pressure positions using
`BMx` and `BMy` (and `BMz`), constructed in operator_divergence.jl.
"""
function get_velocity(V, t, setup)
    (; grid, operators) = setup
    (; Npx, Npy, indu, indv) = grid
    (; Au_ux, Av_vy, Bup, Bvp) = operators

    # Evaluate boundary conditions at current time
    bc_vectors = get_bc_vectors(setup, t)
    (; yAu_ux, yAv_vy) = bc_vectors

    uh = @view V[indu]
    vh = @view V[indv]

    up = reshape(Bup * (Au_ux * uh + yAu_ux), Npx, Npy)
    vp = reshape(Bvp * (Av_vy * vh + yAv_vy), Npx, Npy)

    up, vp
end

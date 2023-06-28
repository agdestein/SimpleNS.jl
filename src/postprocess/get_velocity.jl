"""
    get_velocity(V, setup)

Get velocity values at pressure points. Interpolate velocities to pressure positions using
`BMx` and `BMy` (and `BMz`), constructed in operator_divergence.jl.
"""
function get_velocity(V, setup)
    (; grid, operators) = setup
    (; Npx, Npy, Nx, Ny) = grid
    (; Au_ux, Av_vy, Bup, Bvp) = operators

    u, v = eachslice(reshape(V, Nx, Ny, 2); dims = 3)

    u = reshape(u, :)
    v = reshape(v, :)

    # u = V[indu]
    # v = V[indv]

    up = reshape(Bup * (Au_ux * u), Npx, Npy)
    vp = reshape(Bvp * (Av_vy * v), Npx, Npy)

    up, vp
end

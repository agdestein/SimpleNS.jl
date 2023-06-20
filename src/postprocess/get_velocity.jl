"""
    get_velocity(V, setup)

Get velocity values at pressure points. Interpolate velocities to pressure positions using
`BMx` and `BMy` (and `BMz`), constructed in operator_divergence.jl.
"""
function get_velocity(V, setup)
    (; grid, operators) = setup
    (; Npx, Npy, indu, indv) = grid
    (; Au_ux, Av_vy, Bup, Bvp) = operators

    uh = @view V[indu]
    vh = @view V[indv]

    up = reshape(Bup * (Au_ux * uh), Npx, Npy)
    vp = reshape(Bvp * (Av_vy * vh), Npx, Npy)

    up, vp
end

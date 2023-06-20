"""
    create_body_force(fu, fv, grid)

Steady body force `f(x, y) = [fu(x, y), fv(x, y)]`. 
"""
function create_body_force(fu, fv, grid)
    (; NV, indu, indv, xu, yu, xv, yv) = grid
    F = zeros(eltype(xu), NV)
    F[indu] .= reshape(fu.(xu, yu), :)
    F[indv] .= reshape(fv.(xv, yv), :)
    F
end

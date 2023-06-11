"""
    bodyforce(force, V, t, setup; get_jacobian = false)

Compute body force `F` in momentum equations at velocity points.
"""
function bodyforce end

bodyforce(force::SteadyBodyForce, t, setup) = force.F

function bodyforce(force::UnsteadyBodyForce, t, setup)
    (; indu, indv, xu, xv, yu, yv) = setup.grid
    Fx = force.fu.(xu, yu, t)
    Fy = force.fv.(xv, yv, t)
    [Fx; Fy]
end

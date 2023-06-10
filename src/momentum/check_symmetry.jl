"""
    check_symmetry(V, t, setup, ϵ = 1e-14)

Check symmetry of convection operator.

`flag = 0`: no symmetry error
`flag = 1`: symmetry error
"""
function check_symmetry(V, t, setup, ϵ = 1e-14; bc_vectors)
    (; grid, operators, boundary_conditions) = setup
    (; indu, indv) = grid
    (; Cux, Cuy, Cvx, Cvy) = operators
    (; Au_ux, Au_uy, Av_vx, Av_vy) = operators
    (; Iu_ux, Iv_uy, Iu_vx, Iv_vy) = operators
    (; yIu_ux, yIv_uy, yIu_vx, yIv_vy) = bc_vectors

    uₕ = @view V[indu]
    vₕ = @view V[indv]

    Cu =
        Cux * spdiagm(Iu_ux * uₕ + yIu_ux) * Au_ux +
        Cuy * spdiagm(Iv_uy * vₕ + yIv_uy) * Au_uy
    Cv =
        Cvx * spdiagm(Iu_vx * uₕ + yIu_vx) * Av_vx +
        Cvy * spdiagm(Iv_vy * vₕ + yIv_vy) * Av_vy

    error_u = maximum(abs.(Cu + Cu'))
    error_v = maximum(abs.(Cv + Cv'))

    symmetry_error = max(error_u, error_v)

    flag = 0
    if symmetry_error > ϵ
        if boundary_conditions.u.x[1] != :pressure &&
           boundary_conditions.u.x[2] != :pressure
            flag = 1
        end
        if boundary_conditions.v.y[1] != :pressure &&
           boundary_conditions.v.y[2] != :pressure
            flag = 1
        end
    end

    flag, symmetry_error
end

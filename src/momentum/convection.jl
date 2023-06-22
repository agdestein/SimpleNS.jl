"""
    convection(V, ϕ, setup)

Evaluate convective terms. The convected quantity is `ϕ` (usually `ϕ = V`).
"""
function convection(V, Φ, setup)
    (; grid, operators) = setup
    (; Cux, Cuy, Cvx, Cvy) = operators
    (; Au_ux, Au_uy, Av_vx, Av_vy) = operators
    (; Iu_ux, Iv_uy, Iu_vx, Iv_vy) = operators
    # (; indu, indv) = grid
    (; Nx, Ny) = grid

    u, v = eachslice(reshape(V, Nx, Ny, 2); dims = 3)
    ϕ, ψ = eachslice(reshape(Φ, Nx, Ny, 2); dims = 3)

    u = reshape(u, :)
    v = reshape(v, :)
    ϕ = reshape(ϕ, :)
    ψ = reshape(ψ, :)

    # u = V[indu]
    # v = V[indv]

    # ϕ = Φ[indu]
    # ψ = Φ[indv]

    # Convection components
    u_ux = Au_ux * u                         # u at ux
    ϕ_ux = Iu_ux * ϕ                         # ϕ at ux
    ∂uϕ∂x = Cux * (u_ux .* ϕ_ux)

    u_uy = Au_uy * u                         # u at uy
    ψ_uy = Iv_uy * ψ                         # ψ at uy
    ∂uψ∂y = Cuy * (u_uy .* ψ_uy)

    v_vx = Av_vx * v                         # v at vx
    ϕ_vx = Iu_vx * ϕ                         # ϕ at vx
    ∂vϕ∂x = Cvx * (v_vx .* ϕ_vx)

    v_vy = Av_vy * v                         # v at vy
    ψ_vy = Iv_vy * ψ                         # ψ at vy
    ∂vψ∂y = Cvy * (v_vy .* ψ_vy)

    cu = ∂uϕ∂x + ∂uψ∂y
    cv = ∂vϕ∂x + ∂vψ∂y

    [cu; cv]
end

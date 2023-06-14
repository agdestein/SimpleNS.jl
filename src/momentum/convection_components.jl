"""
    convection_components(
        V, ϕ, setup;
        get_jacobian = false,
        newton_factor = false,
    )

Compute convection components.
"""
function convection_components(
    V,
    Φ,
    setup;
    get_jacobian = false,
    newton_factor = false,
)
    (; grid, operators) = setup
    (; Cux, Cuy, Cvx, Cvy) = operators
    (; Au_ux, Au_uy, Av_vx, Av_vy) = operators
    (; Iu_ux, Iv_uy, Iu_vx, Iv_vy) = operators
    (; indu, indv) = grid

    u = @view V[indu]
    v = @view V[indv]

    ϕ = @view Φ[indu]
    ψ = @view Φ[indv]

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

    c = [cu; cv]

    if get_jacobian
        ## Convective terms, u-component
        C1 = Cux * Diagonal(ϕ_ux)
        C2 = Cux * Diagonal(u_ux) * newton_factor
        Conv_ux_11 = C1 * Au_ux .+ C2 * Iu_ux

        C1 = Cuy * Diagonal(ψ_uy)
        C2 = Cuy * Diagonal(u_uy) * newton_factor
        Conv_uy_11 = C1 * Au_uy
        Conv_uy_12 = C2 * Iv_uy

        ## Convective terms, v-component
        C1 = Cvx * Diagonal(ϕ_vx)
        C2 = Cvx * Diagonal(v_vx) * newton_factor
        Conv_vx_21 = C2 * Iu_vx
        Conv_vx_22 = C1 * Av_vx

        C1 = Cvy * Diagonal(ψ_vy)
        C2 = Cvy * Diagonal(v_vy) * newton_factor
        Conv_vy_22 = C1 * Av_vy .+ C2 * Iv_vy

        ∇c = [
            (Conv_ux_11+Conv_uy_11) Conv_uy_12
            Conv_vx_21 (Conv_vx_22+Conv_vy_22)
        ]
    else
        ∇c = nothing
    end

    c, ∇c
end

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
    ϕ,
    setup;
    get_jacobian = false,
    newton_factor = false,
)
    (; grid, operators) = setup

    (; Cux, Cuy, Cvx, Cvy) = operators
    (; Au_ux, Au_uy, Av_vx, Av_vy) = operators
    (; Iu_ux, Iv_uy, Iu_vx, Iv_vy) = operators

    (; indu, indv) = grid

    uₕ = @view V[indu]
    vₕ = @view V[indv]

    ϕu = @view ϕ[indu]
    ϕv = @view ϕ[indv]

    # Convection components
    u_ux = Au_ux * uₕ                         # u at ux
    ū_ux = Iu_ux * ϕu                         # ū at ux
    ∂uū∂x = Cux * (u_ux .* ū_ux)

    u_uy = Au_uy * uₕ                         # u at uy
    v̄_uy = Iv_uy * ϕv                         # ū at uy
    ∂uv̄∂y = Cuy * (u_uy .* v̄_uy)

    v_vx = Av_vx * vₕ                         # v at vx
    ū_vx = Iu_vx * ϕu                         # ū at vx
    ∂vū∂x = Cvx * (v_vx .* ū_vx)

    v_vy = Av_vy * vₕ                         # v at vy
    v̄_vy = Iv_vy * ϕv                         # ū at vy
    ∂vv̄∂y = Cvy * (v_vy .* v̄_vy)

    cu = ∂uū∂x + ∂uv̄∂y
    cv = ∂vū∂x + ∂vv̄∂y

    c = [cu; cv]

    if get_jacobian
        ## Convective terms, u-component
        C1 = Cux * Diagonal(ū_ux)
        C2 = Cux * Diagonal(u_ux) * newton_factor
        Conv_ux_11 = C1 * Au_ux .+ C2 * Iu_ux

        C1 = Cuy * Diagonal(v̄_uy)
        C2 = Cuy * Diagonal(u_uy) * newton_factor
        Conv_uy_11 = C1 * Au_uy
        Conv_uy_12 = C2 * Iv_uy

        ## Convective terms, v-component
        C1 = Cvx * Diagonal(ū_vx)
        C2 = Cvx * Diagonal(v_vx) * newton_factor
        Conv_vx_21 = C2 * Iu_vx
        Conv_vx_22 = C1 * Av_vx

        C1 = Cvy * Diagonal(v̄_vy)
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

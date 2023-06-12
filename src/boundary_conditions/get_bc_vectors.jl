"""
    get_bc_vectors(setup, t)

Get boundary condition vectors.
"""
function get_bc_vectors(setup::Setup{T,2}, t) where {T}
    (; grid, operators, boundary_conditions, viscosity_model) = setup

    (; Nux_in, Nvy_in, Np, Npx, Npy) = grid
    (; xin, yin, x, y, hx, hy, xp, yp) = grid
    (; Ω, indu, indv) = grid

    (; Dux, Duy, Dvx, Dvy) = operators
    (; Au_ux_bc, Au_uy_bc, Av_vx_bc, Av_vy_bc) = operators
    (; Su_ux_bc, Su_uy_bc, Sv_vx_bc, Sv_vy_bc) = operators
    (; Iu_ux_bc, Iv_uy_bc_lr, Iv_uy_bc_lu) = operators
    (; Iu_vx_bc_lr, Iu_vx_bc_lu, Iv_vy_bc) = operators
    (; Mx_bc, My_bc) = operators
    (; Aν_vy_bc) = operators
    (; Cux_k_bc, Cuy_k_bc, Cvx_k_bc, Cvy_k_bc, Auy_k_bc, Avx_k_bc) = operators
    (; Su_vx_bc_lr, Su_vx_bc_lu, Sv_uy_bc_lr, Sv_uy_bc_lu) = operators

    (; u_bc, v_bc, dudt_bc, dvdt_bc) = boundary_conditions
    (; p_bc) = boundary_conditions

    (; Re) = viscosity_model

    # TODO: Split function into allocating part (constructor?) and mutating `update!`

    ## Get BC values
    uLo = u_bc.(x, y[1], t)
    uUp = u_bc.(x, y[end], t)

    uLo_i = u_bc.(xin, y[1], t)
    uUp_i = u_bc.(xin, y[end], t)
    uLe_i = u_bc.(x[1], yp, t)
    uRi_i = u_bc.(x[end], yp, t)

    vLe = v_bc.(x[1], y, t)
    vRi = v_bc.(x[end], y, t)

    vLo_i = v_bc.(xp, y[1], t)
    vUp_i = v_bc.(xp, y[end], t)
    vLe_i = v_bc.(x[1], yin, t)
    vRi_i = v_bc.(x[end], yin, t)

    ## Boundary conditions for divergence

    # Mx
    ybc = uLe_i ⊗ Mx_bc.ybc1 + uRi_i ⊗ Mx_bc.ybc2
    yMx = Mx_bc.Bbc * ybc

    # My
    ybc = My_bc.ybc1 ⊗ vLo_i + My_bc.ybc2 ⊗ vUp_i
    yMy = My_bc.Bbc * ybc

    yM = yMx + yMy

    # Time derivative of divergence
    ydM = zeros(Np)

    ## Boundary conditions for pressure

    # Left and right side
    y1D_le = zeros(Nux_in)
    y1D_ri = zeros(Nux_in)
    boundary_conditions.u.x[1] == :pressure && (y1D_le[1] = -1)
    boundary_conditions.u.x[2] == :pressure && (y1D_ri[end] = 1)
    y_px = (hy .* p_bc.x[1]) ⊗ y1D_le + (hy .* p_bc.x[2]) ⊗ y1D_ri

    # Lower and upper side
    y1D_lo = zeros(Nvy_in)
    y1D_up = zeros(Nvy_in)
    boundary_conditions.v.y[1] == :pressure && (y1D_lo[1] = -1)
    boundary_conditions.v.y[2] == :pressure && (y1D_up[end] = 1)
    y_py = y1D_lo ⊗ (hx .* p_bc.y[1]) + y1D_up ⊗ (hx .* p_bc.y[2])

    y_p = [y_px; y_py]

    ## Boundary conditions for averaging
    # Au_ux
    ybc = uLe_i ⊗ Au_ux_bc.ybc1 + uRi_i ⊗ Au_ux_bc.ybc2
    yAu_ux = Au_ux_bc.Bbc * ybc

    # Au_uy
    ybc = Au_uy_bc.ybc1 ⊗ uLo_i + Au_uy_bc.ybc2 ⊗ uUp_i
    yAu_uy = Au_uy_bc.Bbc * ybc

    # Av_vx
    ybc = vLe_i ⊗ Av_vx_bc.ybc1 + vRi_i ⊗ Av_vx_bc.ybc2
    yAv_vx = Av_vx_bc.Bbc * ybc

    # Av_vy
    ybc = Av_vy_bc.ybc1 ⊗ vLo_i + Av_vy_bc.ybc2 ⊗ vUp_i
    yAv_vy = Av_vy_bc.Bbc * ybc


    ## Boundary conditions for diffusion
    # Su_ux
    ybc = uLe_i ⊗ Su_ux_bc.ybc1 + uRi_i ⊗ Su_ux_bc.ybc2
    ySu_ux = Su_ux_bc.Bbc * ybc

    # Su_uy
    ybc = Su_uy_bc.ybc1 ⊗ uLo_i + Su_uy_bc.ybc2 ⊗ uUp_i
    ySu_uy = Su_uy_bc.Bbc * ybc

    Sv_uy = Sv_uy_bc_lr.B2D * Sv_uy_bc_lu.B2D

    # Sv_uy (left/right)
    ybc = vLe ⊗ Sv_uy_bc_lr.ybc1 + vRi ⊗ Sv_uy_bc_lr.ybc2
    ySv_uy_lr = Sv_uy_bc_lr.Bbc * ybc

    # Iv_uy (low/up)
    ybc = Sv_uy_bc_lu.ybc1 ⊗ vLo_i + Sv_uy_bc_lu.ybc2 ⊗ vUp_i
    ySv_uy_lu = Sv_uy_bc_lr.B2D * Sv_uy_bc_lu.Bbc * ybc

    ySv_uy = ySv_uy_lr + ySv_uy_lu

    # Su_vx (low/up)
    ybc = Su_vx_bc_lu.ybc1 ⊗ uLo + Su_vx_bc_lu.ybc2 ⊗ uUp
    ySu_vx_lu = Su_vx_bc_lu.Bbc * ybc

    # Su_vx (left/right)
    ybc = uLe_i ⊗ Su_vx_bc_lr.ybc1 + uRi_i ⊗ Su_vx_bc_lr.ybc2
    ySu_vx_lr = Su_vx_bc_lu.B2D * Su_vx_bc_lr.Bbc * ybc
    ySu_vx = ySu_vx_lr + ySu_vx_lu

    # Sv_vx
    ybc = vLe_i ⊗ Sv_vx_bc.ybc1 + vRi_i ⊗ Sv_vx_bc.ybc2
    ySv_vx = Sv_vx_bc.Bbc * ybc

    # Sv_vy
    ybc = Sv_vy_bc.ybc1 ⊗ vLo_i + Sv_vy_bc.ybc2 ⊗ vUp_i
    ySv_vy = Sv_vy_bc.Bbc * ybc

    yDiffu = 1 / Re * (Dux * ySu_ux + Duy * ySu_uy)
    yDiffv = 1 / Re * (Dvx * ySv_vx + Dvy * ySv_vy)
    yDiff = [yDiffu; yDiffv]

    Ωu⁻¹ = 1 ./ Ω[indu]
    Ωv⁻¹ = 1 ./ Ω[indv]

    # Diffusive terms in finite-difference setting, without viscosity
    yDiffu_f = Ωu⁻¹ .* (Dux * ySu_ux + Duy * ySu_uy)
    yDiffv_f = Ωv⁻¹ .* (Dvx * ySv_vx + Dvy * ySv_vy)

    ## Boundary conditions for interpolation

    # Iu_ux
    ybc = uLe_i ⊗ Iu_ux_bc.ybc1 + uRi_i ⊗ Iu_ux_bc.ybc2
    yIu_ux = Iu_ux_bc.Bbc * ybc

    # Iv_uy (left/right)
    ybc = vLe ⊗ Iv_uy_bc_lr.ybc1 + vRi ⊗ Iv_uy_bc_lr.ybc2
    yIv_uy_lr = Iv_uy_bc_lr.Bbc * ybc

    # Iv_uy (low/up)
    ybc = Iv_uy_bc_lu.ybc1 ⊗ vLo_i + Iv_uy_bc_lu.ybc2 ⊗ vUp_i
    yIv_uy_lu = Iv_uy_bc_lr.B2D * Iv_uy_bc_lu.Bbc * ybc
    yIv_uy = yIv_uy_lr + yIv_uy_lu

    # Iu_vx (low/up)
    ybc = Iu_vx_bc_lu.ybc1 ⊗ uLo + Iu_vx_bc_lu.ybc2 ⊗ uUp
    yIu_vx_lu = Iu_vx_bc_lu.Bbc * ybc

    # Iu_vx (left/right)
    ybc = uLe_i ⊗ Iu_vx_bc_lr.ybc1 + uRi_i ⊗ Iu_vx_bc_lr.ybc2
    yIu_vx_lr = Iu_vx_bc_lu.B2D * Iu_vx_bc_lr.Bbc * ybc
    yIu_vx = yIu_vx_lr + yIu_vx_lu

    # Iv_vy
    ybc = Iv_vy_bc.ybc1 ⊗ vLo_i + Iv_vy_bc.ybc2 ⊗ vUp_i
    yIv_vy = Iv_vy_bc.Bbc * ybc

    ## Group BC vectors
    bc_vectors = (;
        yM,
        ydM,
        y_p,
        yAu_ux,
        yAu_uy,
        yAv_vx,
        yAv_vy,
        yDiff,
        yIu_ux,
        yIv_uy,
        yIu_vx,
        yIv_vy,
    )

    # Use values directly (see diffusion.jl)
    bc_vectors = (;
        bc_vectors...,
        ySu_ux,
        ySu_uy,
        ySu_vx,
        ySv_vx,
        ySv_vy,
        ySv_uy,
        yDiffu_f,
        yDiffv_f,
    )

    bc_vectors
end

"""
    get_bc_vectors(setup, t)

Get boundary condition vectors.
"""
function get_bc_vectors end

# 2D version
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
    (; p_bc, bc_unsteady) = boundary_conditions

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

    if bc_unsteady
        dudtLe_i = dudt_bc.(x[1], yp, t)
        dudtRi_i = dudt_bc.(x[end], yp, t)
        dvdtLo_i = dvdt_bc.(xp, y[1], t)
        dvdtUp_i = dvdt_bc.(xp, y[end], t)
    end

    ## Boundary conditions for divergence

    # Mx
    ybc = uLe_i ⊗ Mx_bc.ybc1 + uRi_i ⊗ Mx_bc.ybc2
    yMx = Mx_bc.Bbc * ybc

    # My
    ybc = My_bc.ybc1 ⊗ vLo_i + My_bc.ybc2 ⊗ vUp_i
    yMy = My_bc.Bbc * ybc

    yM = yMx + yMy

    # Time derivative of divergence
    if bc_unsteady
        ybc = dudtLe_i ⊗ Mx_bc.ybc1 + dudtRi_i ⊗ Mx_bc.ybc2
        ydMx = Mx_bc.Bbc * ybc

        # My
        ybc = My_bc.ybc1 ⊗ dvdtLo_i + My_bc.ybc2 ⊗ dvdtUp_i
        ydMy = My_bc.Bbc * ybc

        ydM = ydMx + ydMy
    else
        ydM = zeros(Np)
    end

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

    # Use values directly (see diffusion.jl and strain_tensor.jl)
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

# 3D version
function get_bc_vectors(setup::Setup{T,3}, t) where {T}
    (; grid, operators, boundary_conditions, viscosity_model) = setup

    (; Re) = viscosity_model

    (; u_bc, v_bc, w_bc, dudt_bc, dvdt_bc, dwdt_bc) = boundary_conditions
    (; p_bc, bc_unsteady) = boundary_conditions

    (; Ω, indu, indv, indw) = grid
    (; Nz, Np, Npx, Npy, Npz) = grid
    (; Nux_in, Nux_b, Nux_t, Nuy_in, Nuy_b, Nuy_t, Nuz_in, Nuz_b, Nuz_t) = grid
    (; Nvx_in, Nvx_b, Nvx_t, Nvy_in, Nvy_b, Nvy_t, Nvz_in, Nuz_b, Nvz_t) = grid
    (; Nwx_in, Nwx_b, Nwx_t, Nwy_in, Nwy_b, Nwy_t, Nwz_in, Nwz_b, Nwz_t) = grid
    (; xin, yin, zin, x, y, z, hx, hy, hz, xp, yp, zp) = grid

    (; Dux, Duy, Duz, Dvx, Dvy, Dvz, Dwx, Dwy, Dwz) = operators

    (; Au_ux_bc, Au_uy_bc, Au_uz_bc) = operators
    (; Av_vx_bc, Av_vy_bc, Av_vz_bc) = operators
    (; Aw_wx_bc, Aw_wy_bc, Aw_wz_bc) = operators

    (; Su_ux_bc, Su_uy_bc, Su_uz_bc) = operators
    (; Sv_vx_bc, Sv_vy_bc, Sv_vz_bc) = operators
    (; Sw_wx_bc, Sw_wy_bc, Sw_wz_bc) = operators

    (; Mx_bc, My_bc, Mz_bc) = operators
    (; Aν_vy_bc) = operators

    (; Iu_ux_bc, Iv_vy_bc, Iw_wz_bc) = operators
    (; Iv_uy_bc_lr, Iv_uy_bc_lu, Iw_uz_bc_lr, Iw_uz_bc_bf) = operators
    (; Iu_vx_bc_lr, Iu_vx_bc_lu, Iw_vz_bc_lu, Iw_vz_bc_bf) = operators
    (; Iu_wx_bc_lr, Iu_wx_bc_bf, Iv_wy_bc_lu, Iv_wy_bc_bf) = operators

    (; Cux_k_bc, Cuy_k_bc, Cuz_k_bc) = operators
    (; Cvx_k_bc, Cvy_k_bc, Cvz_k_bc) = operators
    (; Cwx_k_bc, Cwy_k_bc, Cwz_k_bc) = operators

    (; Auy_k_bc, Avx_k_bc) = operators
    (; Auz_k_bc, Awx_k_bc) = operators
    (; Awy_k_bc, Avz_k_bc) = operators

    (; Sv_uy_bc_lr, Sv_uy_bc_lu, Sw_uz_bc_lr, Sw_uz_bc_bf) = operators
    (; Su_vx_bc_lr, Su_vx_bc_lu, Sw_vz_bc_lu, Sw_vz_bc_bf) = operators
    (; Su_wx_bc_lr, Su_wx_bc_bf, Sv_wy_bc_lu, Sv_wy_bc_bf) = operators

    # TODO: Split up function into allocating part (constructor?) and mutating `update!`

    ## Get bc values
    uLe_i = reshape(u_bc.(x[1], yp, zp', t), :)
    uRi_i = reshape(u_bc.(x[end], yp, zp', t), :)
    uLo_i = reshape(u_bc.(xin, y[1], zp', t), :)
    uUp_i = reshape(u_bc.(xin, y[end], zp', t), :)
    uLo_i2 = reshape(u_bc.(x, y[1], zp', t), :)
    uUp_i2 = reshape(u_bc.(x, y[end], zp', t), :)
    uBa_i = reshape(u_bc.(xin, yp', z[1], t), :)
    uFr_i = reshape(u_bc.(xin, yp', z[end], t), :)
    uBa_i2 = reshape(u_bc.(x, yp', z[1], t), :)
    uFr_i2 = reshape(u_bc.(x, yp', z[end], t), :)

    vLe_i = reshape(v_bc.(x[1], yin, zp', t), :)
    vRi_i = reshape(v_bc.(x[end], yin, zp', t), :)
    vLe_i2 = reshape(v_bc.(x[1], y, zp', t), :)
    vRi_i2 = reshape(v_bc.(x[end], y, zp', t), :)
    vLo_i = reshape(v_bc.(xp, y[1], zp', t), :)
    vUp_i = reshape(v_bc.(xp, y[end], zp', t), :)
    vBa_i = reshape(v_bc.(xp, yin', z[1], t), :)
    vFr_i = reshape(v_bc.(xp, yin', z[end], t), :)
    vBa_i2 = reshape(v_bc.(xp, y', z[1], t), :)
    vFr_i2 = reshape(v_bc.(xp, y', z[end], t), :)

    wLe_i = reshape(w_bc.(x[1], yp, zin', t), :)
    wRi_i = reshape(w_bc.(x[end], yp, zin', t), :)
    wLe_i2 = reshape(w_bc.(x[1], yp, z', t), :)
    wRi_i2 = reshape(w_bc.(x[end], yp, z', t), :)
    wLo_i = reshape(w_bc.(xp, y[1], zin', t), :)
    wUp_i = reshape(w_bc.(xp, y[end], zin', t), :)
    wLo_i2 = reshape(w_bc.(xp, y[1], z', t), :)
    wUp_i2 = reshape(w_bc.(xp, y[end], z', t), :)
    wBa_i = reshape(w_bc.(xp, yp', z[1], t), :)
    wFr_i = reshape(w_bc.(xp, yp', z[end], t), :)

    if bc_unsteady
        dudtLe_i = reshape(dudt_bc.(x[1], yp, zp', t), :)
        dudtRi_i = reshape(dudt_bc.(x[end], yp, zp', t), :)
        dvdtLo_i = reshape(dvdt_bc.(xp, y[1], zp', t), :)
        dvdtUp_i = reshape(dvdt_bc.(xp, y[end], zp', t), :)
        dwdtBa_i = reshape(dwdt_bc.(xp, yp', z[1], t), :)
        dwdtFr_i = reshape(dwdt_bc.(xp, yp', z[end], t), :)
    end

    ## Boundary conditions for divergence

    # Mx
    ybc = uLe_i ⊗ Mx_bc.ybc1 + uRi_i ⊗ Mx_bc.ybc2
    yMx = Mx_bc.Bbc * ybc

    # My
    ybc = vLo_i ⊗ My_bc.ybc1 + vUp_i ⊗ My_bc.ybc2
    ybc = reshape(permutedims(reshape(ybc, Nvy_b, Nvx_in, Nvz_in), (2, 1, 3)), :)
    yMy = My_bc.Bbc * ybc

    # Mz
    ybc = Mz_bc.ybc1 ⊗ wBa_i + Mz_bc.ybc2 ⊗ wFr_i
    yMz = Mz_bc.Bbc * ybc

    yM = yMx + yMy + yMz

    # Time derivative of divergence
    if bc_unsteady
        ybc = dudtLe_i ⊗ Mx_bc.ybc1 + dudtRi_i ⊗ Mx_bc.ybc2
        ydMx = Mx_bc.Bbc * ybc

        # My - order of kron is not correct, so reshape
        ybc = dvdtLo_i ⊗ My_bc.ybc1 + dvdtUp_i ⊗ My_bc.ybc2
        ybc = reshape(ybc, Nvy_b, Nvx_in, Nvz_in)
        ybc = permutedims(ybc, (2, 1, 3))
        ybc = reshape(ybc, :)
        ydMy = My_bc.Bbc * ybc

        # Mz
        ybc = Mz_bc.ybc1 ⊗ dwdtBa_i + Mz_bc.ybc2 ⊗ dwdtFr_i
        ydMz = Mz_bc.Bbc * ybc

        ydM = ydMx + ydMy + ydMz
    else
        ydM = zeros(Np)
    end

    ## Boundary conditions for pressure

    # Left and right side
    y1D_le = zeros(Nux_in)
    y1D_ri = zeros(Nux_in)
    boundary_conditions.u.x[1] == :pressure && (y1D_le[1] = -1)
    boundary_conditions.u.x[2] == :pressure && (y1D_ri[end] = 1)
    y_px = (p_bc.x[1] .* (hz ⊗ hy)) ⊗ y1D_le + (p_bc.x[2] .* (hz ⊗ hy)) ⊗ y1D_ri

    # Lower and upper side (order of kron not correct, so reshape)
    y1D_lo = zeros(Nvy_in)
    y1D_up = zeros(Nvy_in)
    boundary_conditions.v.y[1] == :pressure && (y1D_lo[1] = -1)
    boundary_conditions.v.y[2] == :pressure && (y1D_up[end] = 1)
    y_py = (p_bc.y[1] .* (hz ⊗ hx)) ⊗ y1D_lo + (p_bc.y[2] .* (hz ⊗ hx)) ⊗ y1D_up
    y_py = reshape(permutedims(reshape(y_py, Nvy_in, Nvx_in, Nvz_in), (2, 1, 3)), :)

    # Back and front side
    y1D_ba = zeros(Nwz_in)
    y1D_fr = zeros(Nwz_in)
    boundary_conditions.w.z[1] == :pressure && (y1D_ba[1] = -1)
    boundary_conditions.w.z[2] == :pressure && (y1D_fr[end] = 1)
    y_pz = y1D_ba ⊗ (p_bc.z[1] .* (hy ⊗ hx)) + y1D_fr ⊗ (p_bc.z[2] .* (hy ⊗ hx))

    y_p = [y_px; y_py; y_pz]

    ## Boundary conditions for averaging
    # Au_ux
    ybc = uLe_i ⊗ Au_ux_bc.ybc1 + uRi_i ⊗ Au_ux_bc.ybc2
    yAu_ux = Au_ux_bc.Bbc * ybc

    # Au_uy (order of kron is not correct, so permute)
    ybc = Au_uy_bc.ybc1 ⊗ uLo_i + Au_uy_bc.ybc2 ⊗ uUp_i
    ybc = reshape(permutedims(reshape(ybc, Nuy_b, Nux_in, Nuz_in), (2, 3, 1)), :)
    yAu_uy = Au_uy_bc.Bbc * ybc

    # Au_uz
    ybc = Au_uz_bc.ybc1 ⊗ uBa_i + Au_uz_bc.ybc2 ⊗ uFr_i
    yAu_uz = Au_uz_bc.Bbc * ybc

    # Av_vx
    ybc = vLe_i ⊗ Av_vx_bc.ybc1 + vRi_i ⊗ Av_vx_bc.ybc2
    yAv_vx = Av_vx_bc.Bbc * ybc

    # Av_vy (order of kron is not correct, so permute)
    ybc = Av_vy_bc.ybc1 ⊗ vLo_i + Av_vy_bc.ybc2 ⊗ vUp_i
    ybc = reshape(permutedims(reshape(ybc, Nvy_b, Nvx_in, Nvz_in), (2, 3, 1)), :)
    yAv_vy = Av_vy_bc.Bbc * ybc

    # Av_vz
    ybc = Av_vz_bc.ybc1 ⊗ vBa_i + Av_vz_bc.ybc2 ⊗ vFr_i
    yAv_vz = Av_vz_bc.Bbc * ybc

    # Aw_wx
    ybc = wLe_i ⊗ Aw_wx_bc.ybc1 + wRi_i ⊗ Aw_wx_bc.ybc2
    yAw_wx = Aw_wx_bc.Bbc * ybc

    # Aw_wy (order of kron is not correct, so permute)
    ybc = wLo_i ⊗ Aw_wy_bc.ybc1 + wUp_i ⊗ Aw_wy_bc.ybc2
    ybc = reshape(permutedims(reshape(ybc, Nwy_b, Nwx_in, Nwz_in), (2, 3, 1)), :)
    yAw_wy = Aw_wy_bc.Bbc * ybc

    # Aw_wz
    ybc = Aw_wz_bc.ybc1 ⊗ wBa_i + Aw_wz_bc.ybc2 ⊗ wFr_i
    yAw_wz = Aw_wz_bc.Bbc * ybc

    ## Bounary conditions for diffusion

    # Su_ux
    ybc = uLe_i ⊗ Su_ux_bc.ybc1 + uRi_i ⊗ Su_ux_bc.ybc2
    ySu_ux = Su_ux_bc.Bbc * ybc

    # Su_uy
    ybc = uLo_i ⊗ Su_uy_bc.ybc1 + uUp_i ⊗ Su_uy_bc.ybc2
    ybc = reshape(permutedims(reshape(ybc, Nuy_b, Nux_in, Nuz_in), (2, 1, 3)), :)
    ySu_uy = Su_uy_bc.Bbc * ybc

    # Su_uz
    ybc = Su_uz_bc.ybc1 ⊗ uBa_i + Su_uz_bc.ybc2 ⊗ uFr_i
    ySu_uz = Su_uz_bc.Bbc * ybc

    # Sv_uy (left/right)
    ybc = vLe_i2 ⊗ Sv_uy_bc_lr.ybc1 + vRi_i2 ⊗ Sv_uy_bc_lr.ybc2
    ySv_uy_lr = Sv_uy_bc_lr.Bbc * ybc

    # Sv_uy (low/up) (order of kron is not correct, so permute)
    ybc = Sv_uy_bc_lu.ybc1 ⊗ vLo_i + Sv_uy_bc_lu.ybc2 ⊗ vUp_i
    Nb = Nuy_in + 1 - Nvy_in
    Nb != 0 && (ybc = reshape(permutedims(reshape(ybc, Nb, Nvx_in, Nvz_in), (2, 1, 3)), :))
    ySv_uy_lu = Sv_uy_bc_lr.B3D * (Sv_uy_bc_lu.Bbc * ybc)

    # Sv_uy
    ySv_uy = ySv_uy_lr + ySv_uy_lu

    # Sw_uz (left/right)
    ybc = wLe_i2 ⊗ Sw_uz_bc_lr.ybc1 + wRi_i2 ⊗ Sw_uz_bc_lr.ybc2
    ySw_uz_lr = Sw_uz_bc_lr.Bbc * ybc

    # Sw_uz (back/front)
    ybc = Sw_uz_bc_bf.ybc1 ⊗ wBa_i + Sw_uz_bc_bf.ybc2 ⊗ wFr_i
    ySw_uz_bf = Sw_uz_bc_lr.B3D * Sw_uz_bc_bf.Bbc * ybc

    # Sw_uz
    ySw_uz = ySw_uz_lr + ySw_uz_bf

    # Sv_vx
    ybc = vLe_i ⊗ Sv_vx_bc.ybc1 + vRi_i ⊗ Sv_vx_bc.ybc2
    ySv_vx = Sv_vx_bc.Bbc * ybc

    # Sv_vy (order of kron is not correct, so permute)
    ybc = vLo_i ⊗ Sv_vy_bc.ybc1 + vUp_i ⊗ Sv_vy_bc.ybc2
    ybc = reshape(permutedims(reshape(ybc, Nvy_b, Nvx_in, Nvz_in), (2, 1, 3)), :)
    ySv_vy = Sv_vy_bc.Bbc * ybc

    # Sv_vz
    ybc = Sv_vz_bc.ybc1 ⊗ vBa_i + Sv_vz_bc.ybc2 ⊗ vFr_i
    ySv_vz = Sv_vz_bc.Bbc * ybc

    # Su_vx (low/up) (order of kron is not correct, so permute)
    ybc = uLo_i2 ⊗ Su_vx_bc_lu.ybc1 + uUp_i2 ⊗ Su_vx_bc_lu.ybc2
    ybc = reshape(permutedims(reshape(ybc, Nuy_b, Nvx_t - 1, Nvz_in), (2, 1, 3)), :)
    ySu_vx_lu = Su_vx_bc_lu.Bbc * ybc

    # Su_vx (left/right)
    ybc = uLe_i ⊗ Su_vx_bc_lr.ybc1 + uRi_i ⊗ Su_vx_bc_lr.ybc2
    ySu_vx_lr = Su_vx_bc_lu.B3D * Su_vx_bc_lr.Bbc * ybc

    # Su_vx
    ySu_vx = ySu_vx_lr + ySu_vx_lu

    # Sw_vz (low/up) (order of kron is not correct, so permute)
    ybc = wLo_i2 ⊗ Sw_vz_bc_lu.ybc1 + wUp_i2 ⊗ Sw_vz_bc_lu.ybc2
    ybc = reshape(permutedims(reshape(ybc, Nwy_b, Nwx_in, Nz + 1), (2, 1, 3)), :)
    ySw_vz_lu = Sw_vz_bc_lu.Bbc * ybc

    # Sw_vz (back/front)
    ybc = Sw_vz_bc_bf.ybc1 ⊗ wBa_i + Sw_vz_bc_bf.ybc2 ⊗ wFr_i
    ySw_vz_bf = Sw_vz_bc_lu.B3D * Sw_vz_bc_bf.Bbc * ybc

    # Sw_vz
    ySw_vz = ySw_vz_lu + ySw_vz_bf

    # Sw_wx
    ybc = wLe_i ⊗ Sw_wx_bc.ybc1 + wRi_i ⊗ Sw_wx_bc.ybc2
    ySw_wx = Sw_wx_bc.Bbc * ybc

    # Sw_wy (order of kron is not correct, so permute)
    ybc = wLo_i ⊗ Sw_wy_bc.ybc1 + wUp_i ⊗ Sw_wy_bc.ybc2
    ybc = reshape(permutedims(reshape(ybc, Nwy_b, Nwx_in, Nwz_in), (2, 1, 3)), :)
    ySw_wy = Sw_wy_bc.Bbc * ybc

    # Sw_wz
    ybc = Sw_wz_bc.ybc1 ⊗ wBa_i + Sw_wz_bc.ybc2 ⊗ wFr_i
    ySw_wz = Sw_wz_bc.Bbc * ybc

    # Su_wx (back/front)
    ybc = Su_wx_bc_bf.ybc1 ⊗ uBa_i2 + Su_wx_bc_bf.ybc2 ⊗ uFr_i2
    ySu_wx_bf = Su_wx_bc_bf.Bbc * ybc

    # Su_wx (left/right)
    ybc = uLe_i ⊗ Su_wx_bc_lr.ybc1 + uRi_i ⊗ Su_wx_bc_lr.ybc2
    ySu_wx_lr = Su_wx_bc_bf.B3D * Su_wx_bc_lr.Bbc * ybc

    # Su_wx
    ySu_wx = ySu_wx_bf + ySu_wx_lr

    # Sv_wy (back/front)
    ybc = Sv_wy_bc_bf.ybc1 ⊗ vBa_i2 + Sv_wy_bc_bf.ybc2 ⊗ vFr_i2
    ySv_wy_bf = Sv_wy_bc_bf.Bbc * ybc

    # Sv_wy (low/up)
    ybc = vLo_i ⊗ Sv_wy_bc_lu.ybc1 + vUp_i ⊗ Sv_wy_bc_lu.ybc2
    Nb = Nwy_in + 1 - Nvy_in
    Nb != 0 && (ybc = reshape(permutedims(reshape(ybc, Nb, Nvx_in, Nvz_in), (2, 1, 3)), :))
    ySv_wy_lu = Sv_wy_bc_bf.B3D * Sv_wy_bc_lu.Bbc * ybc

    # Sv_wy
    ySv_wy = ySv_wy_bf + ySv_wy_lu

    yDiffu = 1 / Re * (Dux * ySu_ux + Duy * ySu_uy + Duz * ySu_uz)
    yDiffv = 1 / Re * (Dvx * ySv_vx + Dvy * ySv_vy + Dvz * ySv_vz)
    yDiffw = 1 / Re * (Dwx * ySw_wx + Dwy * ySw_wy + Dwz * ySw_wz)
    yDiff = [yDiffu; yDiffv; yDiffw]

    Ωu⁻¹ = 1 ./ Ω[indu]
    Ωv⁻¹ = 1 ./ Ω[indv]
    Ωw⁻¹ = 1 ./ Ω[indw]

    # Diffusive terms in finite-difference setting, without viscosity
    yDiffu_f = Ωu⁻¹ .* (Dux * ySu_ux + Duy * ySu_uy + Duz * ySu_uz)
    yDiffv_f = Ωv⁻¹ .* (Dvx * ySv_vx + Dvy * ySv_vy + Dvz * ySv_vz)
    yDiffw_f = Ωw⁻¹ .* (Dwx * ySw_wx + Dwy * ySw_wy + Dwz * ySw_wz)

    ## Boundary conditions for interpolation

    # Iu_ux
    ybc = uLe_i ⊗ Iu_ux_bc.ybc1 + uRi_i ⊗ Iu_ux_bc.ybc2
    yIu_ux = Iu_ux_bc.Bbc * ybc

    # Iv_uy (left/right)
    ybc = vLe_i2 ⊗ Iv_uy_bc_lr.ybc1 + vRi_i2 ⊗ Iv_uy_bc_lr.ybc2
    yIv_uy_lr = Iv_uy_bc_lr.Bbc * ybc

    # Iv_uy (low/up) (order of kron is not correct, so permute)
    ybc = vLo_i ⊗ Iv_uy_bc_lu.ybc1 + vUp_i ⊗ Iv_uy_bc_lu.ybc2
    Nb = Nuy_in + 1 - Nvy_in
    Nb != 0 && (ybc = reshape(permutedims(reshape(ybc, Nb, Nvx_in, Nvz_in), (2, 1, 3)), :))
    yIv_uy_lu = Iv_uy_bc_lr.B3D * Iv_uy_bc_lu.Bbc * ybc

    yIv_uy = yIv_uy_lr + yIv_uy_lu

    # Iw_uz (left/right)
    ybc = wLe_i2 ⊗ Iw_uz_bc_lr.ybc1 + wRi_i2 ⊗ Iw_uz_bc_lr.ybc2
    yIw_uz_lr = Iw_uz_bc_lr.Bbc * ybc

    # Iw_uz (back/front)
    ybc = Iw_uz_bc_bf.ybc1 ⊗ wBa_i + Iw_uz_bc_bf.ybc2 ⊗ wFr_i
    yIw_uz_bf = Iw_uz_bc_lr.B3D * Iw_uz_bc_bf.Bbc * ybc

    # Iw_uz
    yIw_uz = yIw_uz_lr + yIw_uz_bf

    # Iu_vx (low/up) (order of kron is not correct, so permute)
    ybc = uLo_i2 ⊗ Iu_vx_bc_lu.ybc1 + uUp_i2 ⊗ Iu_vx_bc_lu.ybc2
    ybc = reshape(permutedims(reshape(ybc, Nuy_b, Nvx_t - 1, Nvz_in), (2, 1, 3)), :)
    yIu_vx_lu = Iu_vx_bc_lu.Bbc * ybc

    # Iu_vx (left/right)
    ybc = uLe_i ⊗ Iu_vx_bc_lr.ybc1 + uRi_i ⊗ Iu_vx_bc_lr.ybc2
    yIu_vx_lr = Iu_vx_bc_lu.B3D * Iu_vx_bc_lr.Bbc * ybc

    # Iu_vx
    yIu_vx = yIu_vx_lr + yIu_vx_lu

    # Iv_vy (order of kron is not correct)
    ybc = vLo_i ⊗ Iv_vy_bc.ybc1 + vUp_i ⊗ Iv_vy_bc.ybc2
    ybc = reshape(permutedims(reshape(ybc, Nvy_b, Nvx_in, Nvz_in), (2, 1, 3)), :)
    yIv_vy = Iv_vy_bc.Bbc * ybc

    # Iw_vz (low/up) (order of kron is not correct, so permute)
    ybc = wLo_i2 ⊗ Iw_vz_bc_lu.ybc1 + wUp_i2 ⊗ Iw_vz_bc_lu.ybc2
    ybc = reshape(permutedims(reshape(ybc, Nwy_b, Nwx_in, Nz + 1), (2, 1, 3)), :)
    yIw_vz_lu = Iw_vz_bc_lu.Bbc * ybc

    # Iw_vz (back/front)
    ybc = Iw_vz_bc_bf.ybc1 ⊗ wBa_i + Iw_vz_bc_bf.ybc2 ⊗ wFr_i
    yIw_vz_bf = Iw_vz_bc_lu.B3D * Iw_vz_bc_bf.Bbc * ybc

    # Iw_vz
    yIw_vz = yIw_vz_lu + yIw_vz_bf

    # Iu_wx (back/front)
    ybc = Iu_wx_bc_bf.ybc1 ⊗ uBa_i2 + Iu_wx_bc_bf.ybc2 ⊗ uFr_i2
    yIu_wx_bf = Iu_wx_bc_bf.Bbc * ybc

    # Iu_wx (left/right)
    ybc = uLe_i ⊗ Iu_wx_bc_lr.ybc1 + uRi_i ⊗ Iu_wx_bc_lr.ybc2
    yIu_wx_lr = Iu_wx_bc_bf.B3D * Iu_wx_bc_lr.Bbc * ybc

    yIu_wx = yIu_wx_bf + yIu_wx_lr

    # Iv_wy (back/front)
    ybc = Iv_wy_bc_bf.ybc1 ⊗ vBa_i2 + Iv_wy_bc_bf.ybc2 ⊗ vFr_i2
    yIv_wy_bf = Iv_wy_bc_bf.Bbc * ybc

    # Iv_wy (low/up)
    ybc = vLo_i ⊗ Iv_wy_bc_lu.ybc1 + vUp_i ⊗ Iv_wy_bc_lu.ybc2
    Nb = Nwy_in + 1 - Nvy_in
    Nb != 0 && (ybc = reshape(permutedims(reshape(ybc, Nb, Nvx_in, Nvz_in), (2, 1, 3)), :))
    yIv_wy_lu = Iv_wy_bc_bf.B3D * Iv_wy_bc_lu.Bbc * ybc

    # Iv_wy
    yIv_wy = yIv_wy_bf + yIv_wy_lu

    # Iw_wz
    ybc = Iw_wz_bc.ybc1 ⊗ wBa_i + Iw_wz_bc.ybc2 ⊗ wFr_i
    yIw_wz = Iw_wz_bc.Bbc * ybc

    bc_vectors = (;
        yM,
        ydM,
        y_p,
        yAu_ux,
        yAu_uy,
        yAu_uz,
        yAv_vx,
        yAv_vy,
        yAv_vz,
        yAw_wx,
        yAw_wy,
        yAw_wz,
        yDiff,
        ySu_ux,
        ySu_uy,
        ySu_uz,
        ySv_vx,
        ySv_vy,
        ySv_vz,
        ySw_wx,
        ySw_wy,
        ySw_wz,
        ySu_vx,
        ySu_wx,
        ySv_uy,
        ySv_wy,
        ySw_uz,
        ySw_vz,
        yDiffu_f,
        yDiffv_f,
        yDiffw_f,
        yIu_ux,
        yIv_uy,
        yIw_uz,
        yIu_vx,
        yIv_vy,
        yIw_vz,
        yIu_wx,
        yIv_wy,
        yIw_wz,
    )

    bc_vectors
end

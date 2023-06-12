"""
    Operators()

Discrete operators.
"""
Base.@kwdef struct Operators{T}
    Au_ux::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Au_uy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Au_uz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Av_vx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Av_vy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Av_vz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Aw_wx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Aw_wy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Aw_wz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Au_ux_bc::NamedTuple = (;)
    Au_uy_bc::NamedTuple = (;)
    Au_uz_bc::NamedTuple = (;)
    Av_vx_bc::NamedTuple = (;)
    Av_vy_bc::NamedTuple = (;)
    Av_vz_bc::NamedTuple = (;)
    Aw_wx_bc::NamedTuple = (;)
    Aw_wy_bc::NamedTuple = (;)
    Aw_wz_bc::NamedTuple = (;)

    Iu_ux::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iv_uy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iw_uz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iu_vx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iv_vy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iw_vz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iu_wx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iv_wy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iw_wz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iu_ux_bc::NamedTuple = (;)
    Iv_vy_bc::NamedTuple = (;)
    Iw_wz_bc::NamedTuple = (;)
    Iv_uy_bc_lr::NamedTuple = (;)
    Iv_uy_bc_lu::NamedTuple = (;)
    Iu_vx_bc_lr::NamedTuple = (;)
    Iu_vx_bc_lu::NamedTuple = (;)
    Iw_uz_bc_lr::NamedTuple = (;)
    Iw_uz_bc_bf::NamedTuple = (;)
    Iw_vz_bc_lu::NamedTuple = (;)
    Iw_vz_bc_bf::NamedTuple = (;)
    Iu_wx_bc_lr::NamedTuple = (;)
    Iu_wx_bc_bf::NamedTuple = (;)
    Iv_wy_bc_lu::NamedTuple = (;)
    Iv_wy_bc_bf::NamedTuple = (;)

    M::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Mx_bc::NamedTuple = (;)
    My_bc::NamedTuple = (;)
    Mz_bc::NamedTuple = (;)

    G::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Bup::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Bvp::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Bwp::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Cux::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cuy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cuz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cvx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cvy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cvz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cwx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cwy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cwz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Su_ux::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Su_uy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Su_uz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Su_vx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Su_wx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Sv_vx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Sv_vy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Sv_vz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Sv_uy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Sv_wy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Sw_wx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Sw_wy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Sw_wz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Sw_uz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Sw_vz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Su_ux_bc::NamedTuple = (;)
    Su_uy_bc::NamedTuple = (;)
    Su_uz_bc::NamedTuple = (;)
    Sv_vx_bc::NamedTuple = (;)
    Sv_vy_bc::NamedTuple = (;)
    Sv_vz_bc::NamedTuple = (;)
    Sw_wx_bc::NamedTuple = (;)
    Sw_wy_bc::NamedTuple = (;)
    Sw_wz_bc::NamedTuple = (;)

    Su_vx_bc_lr::NamedTuple = (;)
    Su_vx_bc_lu::NamedTuple = (;)
    Su_wx_bc_lr::NamedTuple = (;)
    Su_wx_bc_bf::NamedTuple = (;)
    Sv_uy_bc_lr::NamedTuple = (;)
    Sv_uy_bc_lu::NamedTuple = (;)
    Sv_wy_bc_lu::NamedTuple = (;)
    Sv_wy_bc_bf::NamedTuple = (;)
    Sw_uz_bc_lr::NamedTuple = (;)
    Sw_uz_bc_bf::NamedTuple = (;)
    Sw_vz_bc_lu::NamedTuple = (;)
    Sw_vz_bc_bf::NamedTuple = (;)

    Dux::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Duy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Duz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Dvx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Dvy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Dvz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Dwx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Dwy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Dwz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Diff::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Wu_uy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Wu_uz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Wv_vx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Wv_vz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Ww_wy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Ww_wx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Cux_k::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cuy_k::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cvx_k::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cvy_k::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Cux_k_bc::NamedTuple = (;)
    Cuy_k_bc::NamedTuple = (;)
    Cuz_k_bc::NamedTuple = (;)
    Cvx_k_bc::NamedTuple = (;)
    Cvy_k_bc::NamedTuple = (;)
    Cvz_k_bc::NamedTuple = (;)
    Cwx_k_bc::NamedTuple = (;)
    Cwy_k_bc::NamedTuple = (;)
    Cwz_k_bc::NamedTuple = (;)

    Auy_k::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Auz_k::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Avx_k::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Avz_k::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Awx_k::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Awy_k::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Auy_k_bc::NamedTuple = (;)
    Auz_k_bc::NamedTuple = (;)
    Avx_k_bc::NamedTuple = (;)
    Avz_k_bc::NamedTuple = (;)
    Awx_k_bc::NamedTuple = (;)
    Awy_k_bc::NamedTuple = (;)

    A::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Aν_ux::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Aν_uy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Aν_uz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Aν_vx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Aν_vy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Aν_vz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Aν_wx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Aν_wy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Aν_wz::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Aν_ux_bc::NamedTuple = (;)
    Aν_vy_bc::NamedTuple = (;)
    Aν_wz_bc::NamedTuple = (;)
    Aν_uy_bc_lr::NamedTuple = (;)
    Aν_uy_bc_lu::NamedTuple = (;)
    Aν_uz_bc_lr::NamedTuple = (;)
    Aν_uz_bc_bf::NamedTuple = (;)
    Aν_vx_bc_lr::NamedTuple = (;)
    Aν_vx_bc_lu::NamedTuple = (;)
    Aν_vz_bc_lu::NamedTuple = (;)
    Aν_vz_bc_bf::NamedTuple = (;)
    Aν_wx_bc_lr::NamedTuple = (;)
    Aν_wx_bc_bf::NamedTuple = (;)
    Aν_wy_bc_lu::NamedTuple = (;)
    Aν_wy_bc_bf::NamedTuple = (;)

    Au_ux3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Au_uy3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Av_vx3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Av_vy3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Au_ux_bc3::NamedTuple = (;)
    Au_uy_bc3::NamedTuple = (;)
    Av_vx_bc3::NamedTuple = (;)
    Av_vy_bc3::NamedTuple = (;)

    Iu_ux3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iv_uy3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iu_vx3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iv_vy3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Iu_ux_bc3::NamedTuple = (;)
    Iv_uy_bc_lr3::NamedTuple = (;)
    Iv_uy_bc_lu3::NamedTuple = (;)
    Iu_vx_bc_lr3::NamedTuple = (;)
    Iu_vx_bc_lu3::NamedTuple = (;)
    Iv_vy_bc3::NamedTuple = (;)

    Mx3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    My3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Mx_bc3::NamedTuple = (;)
    My_bc3::NamedTuple = (;)

    Cux3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cuy3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cvx3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cvy3::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Su_ux_bc3::NamedTuple = (;)
    Su_uy_bc3::NamedTuple = (;)
    Sv_vx_bc3::NamedTuple = (;)
    Sv_vy_bc3::NamedTuple = (;)

    Diffu_f::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Diffv_f::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Diffw_f::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Diffux_div::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Diffuy_div::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Diffvx_div::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Diffvy_div::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
end

"""
    Operators(grid, boundary_conditions, viscosity_model)

Build operators.
"""
function Operators(grid::Grid{T}, boundary_conditions, viscosity_model) where {T}
    # Averaging operators
    op_ave = operator_averaging(grid, boundary_conditions)

    # Interpolation operators
    op_int = operator_interpolation(grid, boundary_conditions)

    # Divergence (u, v) -> p and gradient p -> (u, v) operator
    op_div = operator_divergence(grid, boundary_conditions)

    # Convection operators on u- and v- centered volumes
    op_con = operator_convection_diffusion(grid, boundary_conditions, viscosity_model)

    # Post-processing
    op_pos = operator_postprocessing(grid, boundary_conditions)

    Operators{T}(;
        op_ave...,
        op_int...,
        op_div...,
        op_con...,
        op_reg...,
        op_vis...,
        op_pos...,
    )
end

"""
    Operators()

Discrete operators.
"""
Base.@kwdef struct Operators{T}
    Au_ux::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Au_uy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Av_vx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Av_vy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Iu_ux::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iv_uy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iu_vx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Iv_vy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    M::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    G::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Bup::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Bvp::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Cux::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cuy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cvx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Cvy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Diff::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    Wu_uy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Wv_vx::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    A::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
end

"""
    Operators(grid, viscosity_model)

Build operators.
"""
function Operators(grid::Grid{T}, viscosity_model) where {T}
    # Averaging operators
    op_ave = operator_averaging(grid)

    # Interpolation operators
    op_int = operator_interpolation(grid)

    # Divergence (u, v) -> p and gradient p -> (u, v) operator
    op_div = operator_divergence(grid)

    # Convection operators on u- and v- centered volumes
    op_con = operator_convection_diffusion(grid, viscosity_model)

    # Post-processing
    op_pos = operator_postprocessing(grid)

    Operators{T}(;
        op_ave...,
        op_int...,
        op_div...,
        op_con...,
        op_pos...,
    )
end

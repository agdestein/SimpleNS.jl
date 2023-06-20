"""
    get_operators(grid)

Build operators.
"""
function get_operators(grid)
    # Averaging operators
    op_ave = operator_averaging(grid)

    # Interpolation operators
    op_int = operator_interpolation(grid)

    # Divergence (u, v) -> p and gradient p -> (u, v) operator
    op_div = operator_divergence(grid)

    # Convection operators on u- and v- centered volumes
    op_con = operator_convection_diffusion(grid)

    # Post-processing
    op_pos = operator_postprocessing(grid)

    (;
        op_ave...,
        op_int...,
        op_div...,
        op_con...,
        op_pos...,
    )
end

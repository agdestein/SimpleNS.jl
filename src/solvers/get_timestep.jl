"""
    get_timestep(stepper, cfl; bc_vectors)

Estimate time step based on eigenvalues of operators, using Gershgorin.
"""
function get_timestep(stepper::TimeStepper{M,T,2}, cfl; bc_vectors) where {M,T}
    (; setup, method, V) = stepper
    (; grid, operators) = setup
    (; NV, indu, indv, Ω⁻¹) = grid
    (; Diff) = operators
    (; Cux, Cuy, Cvx, Cvy) = operators
    (; Au_ux, Au_uy, Av_vx, Av_vy) = operators
    (; Iu_ux, Iu_vx, Iv_uy, Iv_vy) = operators
    (; yIu_ux, yIu_vx, yIv_uy, yIv_vy) = bc_vectors

    uₕ = @view V[indu]
    vₕ = @view V[indv]

    # For explicit methods only
    if isexplicit(method)
        # Convective part
        Cu =
            Cux * spdiagm(Iu_ux * uₕ + yIu_ux) * Au_ux +
            Cuy * spdiagm(Iv_uy * vₕ + yIv_uy) * Au_uy
        Cv =
            Cvx * spdiagm(Iu_vx * uₕ + yIu_vx) * Av_vx +
            Cvy * spdiagm(Iv_vy * vₕ + yIv_vy) * Av_vy
        C = blockdiag(Cu, Cv)
        test = spdiagm(Ω⁻¹) * C
        sum_conv = abs.(test) * ones(NV) - diag(abs.(test)) - diag(test)
        λ_conv = maximum(sum_conv)

        # Based on max. value of stability region (not a very good indication
        # For the methods that do not include the imaginary axis)
        Δt_conv = lambda_conv_max(method) / λ_conv

        ## Diffusive part
        test = Diagonal(Ω⁻¹) * Diff
        sum_diff = abs.(test) * ones(NV) - diag(abs.(test)) - diag(test)
        λ_diff = maximum(sum_diff)

        # Based on max. value of stability region
        Δt_diff = lambda_diff_max(method) / λ_diff

        Δt = cfl * min(Δt_conv, Δt_diff)
    end

    Δt
end

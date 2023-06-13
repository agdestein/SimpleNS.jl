"""
    diffusion!(model, V, setup; get_jacobian = false)

Evaluate diffusive terms `d` and optionally Jacobian `∇d = ∂d/∂V` using viscosity model `model`.
"""
function diffusion(::LaminarModel, V, setup; get_jacobian = false)
    (; Diff) = setup.operators

    d = Diff * V

    if get_jacobian
        ∇d = Diff
    else
        ∇d = nothing
    end

    d, ∇d
end

"""
    diffusion!(model, V, setup; bc_vectors, get_jacobian = false)

Evaluate diffusive terms `d` and optionally Jacobian `∇d = ∂d/∂V` using viscosity model `model`.
"""
function diffusion(::LaminarModel, V, setup; bc_vectors, get_jacobian = false)
    (; Diff) = setup.operators
    (; yDiff) = bc_vectors

    d = Diff * V + yDiff

    if get_jacobian
        ∇d = Diff
    else
        ∇d = nothing
    end

    d, ∇d
end

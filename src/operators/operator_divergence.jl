"""
    operator_divergence(grid)

Construct divergence and gradient operator.
"""
function operator_divergence(grid)
    (; Npx, Npy) = grid
    (; Nux_in, Nux_b, Nux_t, Nuy_in) = grid
    (; Nvx_in, Nvy_in, Nvy_b, Nvy_t) = grid
    (; hx, hy) = grid
    (; Ω⁻¹) = grid

    ## Divergence operator M

    # Note that the divergence matrix M is not square
    mat_hx = Diagonal(hx)
    mat_hy = Diagonal(hy)

    # For fourth order: mat_hx3 is defined in operator_interpolation

    ## Mx
    # Building blocks consisting of diagonal matrices where the diagonal is
    # equal to constant per block (hy(block)) and changing for next block to
    # hy(block+1)
    diag1 = ones(Nux_t - 1)
    M1D = spdiagm(Nux_t - 1, Nux_t, 0 => -diag1, 1 => diag1)

    # We only need derivative at inner pressure points, so we map the resulting
    # boundary matrix (restrict)
    diagpos = 1

    BMx = spdiagm(Npx, Nux_t - 1, diagpos => ones(Npx))
    M1D = BMx * M1D

    # Extension to 2D to be used in post-processing files
    Bup = I(Nuy_in) ⊗ BMx

    # Boundary conditions
    Mx_bc = bc_general(Nux_t, Nux_in, Nux_b, hx[1], hx[end])
    Mx_bc = (; Mx_bc..., Bbc = mat_hy ⊗ (M1D * Mx_bc.Btemp))

    # Extend to 2D
    Mx = mat_hy ⊗ (M1D * Mx_bc.B1D)

    ## My (same as Mx but reversing indices and kron arguments)
    diag1 = ones(Nvy_t - 1)
    M1D = spdiagm(Nvy_t - 1, Nvy_t, 0 => -diag1, 1 => diag1)

    # We only need derivative at inner pressure points, so we map the resulting
    # boundary matrix (restriction)
    diagpos = 1

    BMy = spdiagm(Npy, Nvy_t - 1, diagpos => ones(Npy))
    M1D = BMy * M1D

    # Extension to 2D to be used in post-processing files
    Bvp = BMy ⊗ I(Nvx_in)

    # Boundary conditions
    My_bc = bc_general(Nvy_t, Nvy_in, Nvy_b, hy[1], hy[end])
    My_bc = (; My_bc..., Bbc = (M1D * My_bc.Btemp) ⊗ mat_hx)

    # Extend to 2D
    My = (M1D * My_bc.B1D) ⊗ mat_hx

    ## Resulting divergence matrix
    M = [Mx My]

    ## Gradient operator G

    # Like in the continuous case, grad = -div^T
    # Note that this also holds for outflow boundary conditions, if the stress
    # on the ouflow boundary is properly taken into account in y_p (often this
    # stress will be zero)
    G = -M'

    ## Pressure matrix for pressure correction method;
    # Also used to make initial data divergence free or compute additional poisson solve
    # Note that the matrix for the pressure is constant in time.
    # Only the right hand side vector changes, so the pressure matrix can be set up outside the time-stepping-loop.

    # Laplace = div grad
    A = M * Diagonal(Ω⁻¹) * G

    # Check if all the row sums of the pressure matrix are zero, which
    # should be the case if there are no pressure boundary conditions
    if any(≉(0; atol = 1e-10), sum(A; dims = 2))
        @warn "Pressure matrix: not all rowsums are zero!"
    end

    ## Group operators
    operators = (; M, G, Bup, Bvp, A)

    operators
end

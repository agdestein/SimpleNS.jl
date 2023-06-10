"""
    strain_tensor(V, setup; bc_vectors, get_jacobian = false, get_S_abs = false)

Evaluate rate of strain tensor `S(V)` and its magnitude.
"""
function strain_tensor end

# 2D version
function strain_tensor(
    V,
    setup::Setup{T,2};
    bc_vectors,
    get_jacobian = false,
    get_S_abs = false,
) where {T}
    (; grid, operators, boundary_conditions) = setup
    (; Nx, Ny, Nu, Nv, Np, indu, indv) = grid
    (; Nux_in, Nuy_in, Nvx_in, Nvy_in) = grid
    (; x, y, xp, yp) = grid
    (; Su_ux, Su_uy, Su_vx, Sv_vx, Sv_vy, Sv_uy) = operators
    (; Cux_k, Cuy_k, Cvx_k, Cvy_k, Auy_k, Avx_k) = operators
    (; ySu_ux, ySu_uy, ySu_vx, ySv_vx, ySv_vy, ySv_uy) = bc_vectors
    (; yCux_k, yCuy_k, yCvx_k, yCvy_k, yAuy_k, yAvx_k) = bc_vectors

    uₕ = @view V[indu]
    vₕ = @view V[indv]

    # These four components are approximated by
    S11 = Su_ux * uₕ + ySu_ux
    S12 = 1 / 2 * (Su_uy * uₕ + ySu_uy + Sv_uy * vₕ + ySv_uy)
    S21 = 1 / 2 * (Su_vx * uₕ + ySu_vx + Sv_vx * vₕ + ySv_vx)
    S22 = Sv_vy * vₕ + ySv_vy

    # Note: S11 and S22 at xp, yp locations (pressure locations)
    # S12, S21 at vorticity locations (corners of pressure cells, (x, y))

    # Option 1: get each S11, S12, S21, S22 at 4 locations (ux locations, at uy locations, at vx
    # Locations and at vy locations); this gives 16 S fields. determine S_abs at each of these
    # Locations, giving 4 S_abs fields, that can be used in computing
    # Dux*(S_abs_ux .* (Su_ux*uₕ+ySu_ux)) etc.

    # Option 2: interpolate S11, S12, S21, S22 to pressure locations
    # Determine S_abs at pressure location
    # Then interpolate to ux, uy, vx, vy locations

    # We will use option 2;
    # Within option 2, we can decide to interpolate S12, S21 etc (option 2a), or we can use
    # Directly the operators that map from velocity field to the S locations,
    # As used for example in ke_production (option 2b).

    if get_S_abs
        # Option 2b
        if boundary_conditions.u.x == (:periodic, :periodic)
            # "cut-off" the double points in case of periodic BC
            # For periodic boundary conditions S11(Npx+1, :) = S11(1, :)
            # So S11 has size (Npx+1)*Npy; the last row are "ghost" points equal to the
            # First points. we have S11 at positions ([xp[1] - 1/2*(hx[1]+hx[end]); xp], yp)
            S11_p = reshape(S11, Nux_in + 1, Nuy_in)
            S11_p = S11_p(2:(Nux_in+1), :) # B

            # S12 is defined on the corners: size Nux_in*(Nuy_in+1), positions (xin, y)
            # Get S12 and S21 at all corner points
            S12_temp = zeros(Nx + 1, Ny + 1)
            S12_temp[1:Nx, :] = reshape(S12, Nx, Ny + 1)
            S12_temp[Nx+1, :] = S12_temp[1, :]
        elseif boundary_conditions.u.x[1] == :dirichlet &&
               boundary_conditions.u.x[2] == :pressure
            S11_p = reshape(S11, Nux_in + 1, Nuy_in)
            S11_p = S11_p[1:Nux_in, :] # Cut off last point

            # S12 is defined on the corners: size Nux_in*(Nuy_in+1), positions (xin, y)
            # Get S12 and S21 at all corner points
            S12_temp = zeros(Nx + 1, Ny + 1)
            S12_temp[2:(Nx+1), :] = reshape(S12, Nx, Ny + 1)
            S12_temp[1, :] = S12_temp[2, :] # Copy from x[2] to x[1]; one could do this more accurately in principle by using the BC
        else
            error("BC not implemented in strain_tensor.jl")
        end

        if boundary_conditions.v.y[1] == :periodic &&
           boundary_conditions.v.y[2] == :periodic
            # Similarly, S22(:, Npy+1) = S22(:, 1). positions (xp, [yp;yp[1]])
            S22_p = reshape(S22, Nvx_in, Nvy_in + 1)
            S22_p = S22_p(:, 2:(Nvy_in+1)) # Why not 1:Nvy_in?

            # Similarly S21 is size (Nux_in+1)*Nuy_in, positions (x, yin)
            S21_temp = zeros(Nx + 1, Ny + 1)
            S21_temp[:, 1:Ny] = reshape(S21, Nx + 1, Ny)
            S21_temp[:, Ny+1] = S21_temp[:, 1]
        elseif boundary_conditions.v.y[1] == :pressure &&
               boundary_conditions.v.y[2] == :pressure
            S22_p = reshape(S22, Nvx_in, Nvy_in + 1)
            S22_p = S22_p(:, 2:Nvy_in)

            # This is nicely defined on all corners
            S21_temp = reshape(S21, Nx + 1, Ny + 1)
        else
            error("BC not implemented in strain_tensor.jl")
        end

        # Now interpolate S12 and S21 to pressure points
        # S11 and S22 have already been trimmed down to this grid

        # S21 and S12 should be equal!
        # FIXME: Find interpolation syntax
        error("Interpolation not implemented")
        # S12_p = interp2(y', x, S12_temp, yp', xp)
        # S21_p = interp2(y', x, S21_temp, yp', xp)

        ## Invariants
        q = @. 1 / 2 * (S11_p[:]^2 + S12_p[:]^2 + S21_p[:]^2 + S22_p[:]^2)

        # Absolute value of strain tensor
        # With S as defined above, i.e. 1/2*(grad u + grad u^T)
        # S_abs = sqrt(2*tr(S^2)) = sqrt(4*q)
        S_abs = sqrt(4q)
    else
        # Option 2a
        S11_p = Cux_k * uₕ + yCux_k
        S12_p =
            1 / 2 * (
                Cuy_k * (Auy_k * uₕ + yAuy_k) +
                yCuy_k +
                Cvx_k * (Avx_k * vₕ + yAvx_k) +
                yCvx_k
            )
        S21_p = S12_p
        S22_p = Cvy_k * vₕ + yCvy_k

        S_abs = @. sqrt(2 * S11_p^2 + 2 * S22_p^2 + 2 * S12_p^2 + 2 * S21_p^2)

        # Jacobian of S_abs wrt u and v
        if get_jacobian
            eps = 1e-14
            Sabs_inv = spdiagm(1 ./ (2 .* S_abs .+ eps))
            Jacu =
                Sabs_inv * (4 * spdiagm(S11_p) * Cux_k + 4 * spdiagm(S12_p) * Cuy_k * Auy_k)
            Jacv =
                Sabs_inv * (4 * spdiagm(S12_p) * Cvx_k * Avx_k + 4 * spdiagm(S22_p) * Cvy_k)
        else
            Jacu = spzeros(Np, Nu)
            Jacv = spzeros(Np, Nv)
        end
    end

    S11, S12, S21, S22, S_abs, Jacu, Jacv
end

# 3D version
function strain_tensor(
    V,
    setup::Setup{T,3};
    bc_vectors,
    get_jacobian = false,
    get_S_abs = false,
) where {T}
    error("Not implemented (3D)")
end

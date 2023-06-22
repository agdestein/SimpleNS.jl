"""
Total solution u is written as `u = Bb*ub + Bin*uin`
The boundary conditions can be written as
`Bbc*u = ybc`
Then u can be written entirely in terms of uin and ybc as:
`u = (Bin-Btemp*Bbc*Bin)*uin + Btemp*ybc`,
where
`Btemp = Bb*(Bbc*Bb)^(-1)`
`Bb`, `Bin` and `Bbc` depend on type of BC (Neumann/Dirichlet/periodic)
Val1 and val2 can be scalars or vectors with either the value or the
derivative
(ghost) points on staggered locations (pressure lines)
"""
function bc_general_stag(T, Nt, Nin, Nb)

    # @show Nt Nin Nb

    # Some input checking:
    if Nt != Nin + Nb
        error("Number of inner points plus boundary points is not equal to total points")
    end

    # Boundary conditions
    Bbc = spzeros(T, Nb, Nt)

    if Nb == 0
        # No boundary points, so simply diagonal matrix without boundary contribution
        B1D = I(Nt)
        Btemp = spzeros(T, Nt, 2)
    elseif Nb ∈ (1, 2)
        if Nb == 1
            # One boundary point
            Bb = spzeros(T, Nt, Nb)

            diagpos = -1
            Bbc[1, 1] = -1
            Bbc[1, end] = 1
            Bb[1, 1] = 1

            # Boundary matrices
            Bin = spdiagm(Nt, Nin, diagpos => ones(T, Nin))
        elseif Nb == 2
            # Normal situation, 2 boundary points
            # Boundary matrices
            Bin = spdiagm(Nt, Nin, -1 => ones(T, Nin))
            Bb = spzeros(T, Nt, Nb)
            Bb[1, 1] = 1
            Bb[end, Nb] = 1

            Bbc[1, 1] = -1
            Bbc[1, end-1] = 1
            Bbc[2, 2] = -1
            Bbc[2, end] = 1
        end
        Btemp = Bb / (Bbc * Bb)
        B1D = Bin - Btemp * Bbc * Bin
    else
        error("Nb must be 0, 1, or 2")
    end

    (; B1D, Btemp)
end

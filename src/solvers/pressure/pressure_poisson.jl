"""
    pressure_poisson(solver, f, setup)

Solve the Poisson equation for the pressure with right hand side `f` at time `t`.
For periodic and no-slip BC, the sum of `f` should be zero.
"""
function pressure_poisson end

function pressure_poisson(solver::DirectPressureSolver, f)
    # Assume the Laplace matrix is known (A) and is possibly factorized

    # Use pre-determined decomposition
    solver.A_fact \ f
end

function pressure_poisson(solver::CGPressureSolver, f)
    (; A, abstol, reltol, maxiter) = solver
    cg(A, f; abstol, reltol, maxiter)
end

function pressure_poisson(solver::FourierPressureSolver, f)
    (; Ahat) = solver

    f = reshape(f, size(Ahat))

    # Fourier transform of right hand side
    fhat = fft(f)

    # Solve for coefficients in Fourier space
    phat = @. -fhat / Ahat

    # Transform back
    p = ifft(phat)

    reshape(real.(p), :)
end

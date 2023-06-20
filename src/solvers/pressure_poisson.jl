"""
    pressure_poisson(setup, f)

Solve the Poisson equation for the pressure with right hand side `f` at time `t`.
For periodic and no-slip BC, the sum of `f` should be zero.
"""
function pressure_poisson(setup, f)
    (; Ahat) = setup.operators

    f = reshape(f, size(Ahat))

    # Fourier transform of right hand side
    fhat = fft(f)

    # Solve for coefficients in Fourier space
    phat = @. -fhat / Ahat

    # Transform back
    p = ifft(phat)

    reshape(real.(p), :)
end

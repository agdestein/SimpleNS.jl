"""
    function solve(
        V, p, tlims;
        setup,
        method = RK44(),
        Δt = nothing,
        processors = [],
    )

Solve unsteady problem using `method`.

If `Δt` is a real number, it is rounded such that `(t_end - t_start) / Δt` is
an integer.
If `Δt = nothing`, the time step is chosen every `n_adapt_Δt` iteration with
CFL-number `cfl` .

For methods that are not self-starting, `nstartup` startup iterations are performed with
`method_startup`.

Each `processor` is called after every `processor.nupdate` time step.
"""
function solve(V, p, tlims; setup, method = RK44(), Δt = nothing, processors = [])
    t_start, t_end = tlims
    nstep = round(Int, (t_end - t_start) / Δt)
    Δt = (t_end - t_start) / nstep

    stepper = (; method, setup, V, p, t = t_start, n = 1)

    # Processors for iteration results  
    for ps ∈ processors
        initialize!(ps, stepper)
        process!(ps, stepper)
    end

    for i = 1:nstep
        # Perform a single time step with the time integration method
        # stepper = step(stepper, Δt)
        stepper = step_rk4(stepper, Δt)

        # Process iteration results with each processor
        for ps ∈ processors
            # Only update each `nupdate`-th iteration
            stepper.n % ps.nupdate == 0 && process!(ps, stepper)
        end
    end

    foreach(finalize!, processors)

    (; V, p) = stepper
    V, p
end

function create_data(
    tburn,
    tlims;
    setup,
    M,
    method = RK44(),
    Δt = nothing,
    nsubstep = 10,
    device = identity,
)
    T = eltype(tlims)
    (; N) = setup.grid
    t_start, t_end = tlims
    nstep = round(Int, (t_end - t_start) / Δt)
    Δt = (t_end - t_start) / nstep

    nburn = round(Int, tburn / Δt)

    nsave = nstep ÷ nsubstep
    @assert nsave * nsubstep == nstep

    # Initial conditions (on device)
    K = N ÷ 2
    V, p = random_field(
        setup,
        K;
        A = T(1e8),
        σ = T(30),
        ## σ = 10,
        s = 5,
        device,
    )

    # Filter matrix (on device)
    W = create_top_hat_velocity(N, M)

    # Output arrays (on host)
    Vbar = zeros(T, 2 * M^2, nsave + 1)
    Fbar = zeros(T, 2 * M^2, nsave + 1)
    tbar = LinRange(t_start, t_end, nsave + 1)
    t = T(0)

    stepper = (; method, setup, V, p, t = t_start, n = 1);

    # Do some burn-in steps
    @info "Running burn in simulation for $tburn seconds"
    for j = 1:nburn
        @show t
        stepper = step_rk4(stepper, Δt)
        t += Δt
    end

    t = t_start


    @info "Running DNS simulation for $(t_end - t_start) seconds"
    for i = 1:nsave
        @show t
        (; V, p) = stepper
        Vbar[:, i] = apply_matrix(W, V)

        F = momentum(V, V, p, setup; nopressure = true)
        Fbar[:, i] = apply_matrix(W, F)

        for j = 1:nsubstep
            # Perform a single time step with the time integration method
            # stepper = step(stepper, Δt)
            stepper = step_rk4(stepper, Δt)
            t += Δt
        end
    end

    (; V, p) = stepper

    Vbar[:, end] = apply_matrix(W, V)

    F = momentum(V, V, p, setup; nopressure = true)
    Fbar[:, end] = apply_matrix(W, F)

    (; Vbar, Fbar, tbar, nsubstep)
end

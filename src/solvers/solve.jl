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
function solve(
    V, p, tlims;
    setup,
    method = RK44(),
    Δt = nothing,
    processors = [],
)
    t_start, t_end = tlims
    nstep = round(Int, (t_end - t_start) / Δt)
    Δt = (t_end - t_start) / nstep

    stepper = (;
        method,
        setup,
        V,
        p,
        t = t_start,
        n = 1,
    )

    # Processors for iteration results  
    for ps ∈ processors
        initialize!(ps, stepper)
        process!(ps, stepper)
    end

    for i = 1:nstep
        # Perform a single time step with the time integration method
        stepper = step(stepper, Δt)

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

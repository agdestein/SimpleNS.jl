"""
    plot_velocity(setup, V; kwargs...)

Plot velocity.
"""
function plot_velocity(setup, V; kwargs...)
    (; xp, yp, xlims, ylims) = setup.grid

    # Get velocity at pressure points
    up, vp = get_velocity(V, setup)
    qp = map((u, v) -> √(u^2 + v^2), up, vp)

    # Levels
    μ, σ = mean(qp), std(qp)
    ≈(μ + σ, μ; rtol = 1e-8, atol = 1e-8) && (σ = 1e-4)
    levels = LinRange(μ - 1.5σ, μ + 1.5σ, 10)

    fig = Figure()
    ax = Axis(
        fig[1, 1];
        aspect = DataAspect(),
        title = "Velocity",
        xlabel = "x",
        ylabel = "y",
    )
    limits!(ax, xlims[1], xlims[2], ylims[1], ylims[2])
    cf = contourf!(ax, xp, yp, qp; extendlow = :auto, extendhigh = :auto, levels, kwargs...)
    Colorbar(fig[1, 2], cf)
    # Colorbar(fig[2,1], cf; vertical = false)
    fig
end

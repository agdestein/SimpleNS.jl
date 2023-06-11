"""
    plot_force(setup, t; kwargs...)

Plot body force.
"""
function plot_force(setup, t; kwargs...)
    (; grid, force) = setup
    (; xp, yp, xlims, ylims) = grid
    (; xu, yu, zu, indu, xlims, ylims) = grid
    (; F) = force

    Fp = reshape(F[indu], size(xu))
    # TODO: norm of F instead of Fu
    xp, yp = xu[:, 1], yu[1, :]

    # Levels
    μ, σ = mean(Fp), std(Fp)
    ≈(μ + σ, μ; rtol = 1e-8, atol = 1e-8) && (σ = 1e-4)
    levels = LinRange(μ - 1.5σ, μ + 1.5σ, 10)

    fig = Figure()
    ax = Axis(fig[1, 1]; aspect = DataAspect(), title = "Force", xlabel = "x", ylabel = "y")
    limits!(ax, xlims[1], xlims[2], ylims[1], ylims[2])
    cf = contourf!(ax, xp, yp, Fp; extendlow = :auto, extendhigh = :auto, levels, kwargs...)
    Colorbar(fig[1, 2], cf)
    fig
end

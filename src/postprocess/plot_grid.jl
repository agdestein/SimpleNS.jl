"""
    plot_grid(x, y)
    plot_grid(grid)

Plot nonuniform Cartesian grid.
"""
function plot_grid end

plot_grid(g::Grid) = plot_grid(g.x, g.y)

plot_grid(x, y) = wireframe(
    x,
    y,
    zeros(length(x), length(y));
    axis = (; aspect = DataAspect(), xlabel = "x", ylabel = "y"),
)

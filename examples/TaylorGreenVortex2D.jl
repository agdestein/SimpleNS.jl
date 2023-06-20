# Little LSP hack to get function signatures, go    #src
# to definition etc.                                #src
if isdefined(@__MODULE__, :LanguageServer)          #src
    include("../src/IncompressibleNavierStokes.jl") #src
    using .IncompressibleNavierStokes               #src
end                                                 #src

# # Taylor-Green vortex - 2D
#
# In this example we consider the Taylor-Green vortex.

# We start by loading packages.
# A [Makie](https://github.com/JuliaPlots/Makie.jl) plotting backend is needed
# for plotting. `GLMakie` creates an interactive window (useful for real-time
# plotting), but does not work when building this example on GitHub.
# `CairoMakie` makes high-quality static vector-graphics plots.

#md using CairoMakie
using GLMakie #!md
using IncompressibleNavierStokes

# Case name for saving results
name = "TaylorGreenVortex2D"

# Viscosity model
viscosity_model = LaminarModel(; Re = 2000.0)

# A 2D grid is a Cartesian product of two vectors
n = 100
x = LinRange(0, 2π, n + 1)
y = LinRange(0, 2π, n + 1)
plot_grid(x, y)

# Build setup and assemble operators
setup = get_setup(x, y; viscosity_model);

# Time interval
t_start, t_end = tlims = (0.0, 10.0)

# Initial conditions
initial_velocity_u(x, y) = -sin(x)cos(y)
initial_velocity_v(x, y) = cos(x)sin(y)
initial_pressure(x, y) = 1 / 4 * (cos(2x) + cos(2y))
V₀, p₀ = create_initial_conditions(
    setup;
    initial_velocity_u,
    initial_velocity_v,
    initial_pressure,
);

# Iteration processors
logger = Logger()
observer = StateObserver(1, V₀, p₀, t_start)
writer = VTKWriter(; nupdate = 10, dir = "output/$name", filename = "solution")
## processors = [logger, observer, writer]
processors = [logger, observer]

# Real time plot
real_time_plot(observer, setup)

# Solve unsteady problem
V, p = solve(V₀, p₀, tlims; setup, Δt = 0.01, processors)
#md current_figure()

# ## Post-process
#
# We may visualize or export the computed fields `(V, p)`

# Export to VTK
save_vtk(V, p, setup, "output/solution")

# Plot pressure
plot_pressure(setup, p)

# Plot velocity
plot_velocity(setup, V)

# Plot vorticity
plot_vorticity(setup, V)

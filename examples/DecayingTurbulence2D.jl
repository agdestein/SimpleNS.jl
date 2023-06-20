# Little LSP hack to get function signatures, go    #src
# to definition etc.                                #src
if isdefined(@__MODULE__, :LanguageServer)          #src
    include("../src/IncompressibleNavierStokes.jl") #src
    using .IncompressibleNavierStokes               #src
end                                                 #src

# # Decaying Homogeneous Isotropic Turbulence - 2D
#
# In this example we consider decaying homogeneous isotropic turbulence,
# similar to the cases considered in [Kochkov2021](@cite) and
# [Kurz2022](@cite). The initial velocity field is created randomly, but with a
# specific energy spectrum. Due to viscous dissipation, the turbulent features
# eventually group to form larger visible eddies.

# We start by loading packages.
# A [Makie](https://github.com/JuliaPlots/Makie.jl) plotting backend is needed
# for plotting. `GLMakie` creates an interactive window (useful for real-time
# plotting), but does not work when building this example on GitHub.
# `CairoMakie` makes high-quality static vector-graphics plots.

using CUDA
using FFTW
#md using CairoMakie
using GLMakie #!md
using IncompressibleNavierStokes
using LaTeXStrings
using LinearAlgebra

# Case name for saving results
name = "DecayingTurbulence2D"

# A 2D grid is a Cartesian product of two vectors
Nx = 256
Ny = 256
xlims = (0.0, 1.0)
ylims = (0.0, 1.0)

# Viscosity
Re = 1e4
viscosity = 1 / Re

# Build setup and assemble operators
setup = get_setup(Nx, Ny, xlims, ylims; viscosity);

# Initial conditions
K = Nx ÷ 2
A = 1e6
σ = 30
## σ = 10
s = 5
function create_spectrum(K)
    a =
        A * [
            1 / sqrt((2π)^2 * 2σ^2) *
            exp(-((i - s)^2 + (j - s)^2) / 2σ^2) *
            exp(-2π * im * rand()) for i = 1:K, j = 1:K
        ]
    [
        a reverse(a; dims = 2)
        reverse(a; dims = 1) reverse(a)
    ]
end
u = real.(ifft(create_spectrum(K)))
v = real.(ifft(create_spectrum(K)))
V = [reshape(u, :); reshape(v, :)]
f = setup.operators.M * V
p = zero(f)

# Make velocity field divergence free
(; Ω) = setup.grid;
(; G, M) = setup.operators;
f = M * V;
Δp = pressure_poisson(setup, f);
V .-= (G * Δp ./  Ω);
p = pressure_additional_solve(setup, V, p);

V₀, p₀ = V, p

# Iteration processors
logger = Logger()
observer = StateObserver(1, V₀, p₀, 0.0)
writer = VTKWriter(; nupdate = 10, dir = "output/$name", filename = "solution")
# processors = [logger, observer, writer]
processors = [logger, writer]
# processors = [logger, observer]
# processors = [logger]

# Real time plot
rtp = real_time_plot(observer, setup; type = heatmap)

# Plot energy history
(; Ωp) = setup.grid
_points = Point2f[]
points = @lift begin
    V, p, t = $(observer.state)
    up, vp = get_velocity(V, t, setup)
    up = reshape(up, :)
    vp = reshape(vp, :)
    E = up' * Diagonal(Ωp) * up + vp' * Diagonal(Ωp) * vp
    push!(_points, Point2f(t, E))
end
ehist = lines(points; axis = (; xlabel = "t", ylabel = "Kinetic energy"))

# Plot energy spectrum
k = 1:(K-1)
kk = reshape([sqrt(kx^2 + ky^2) for kx ∈ k, ky ∈ k], :)
ehat = @lift begin
    V, p, t = $(observer.state)
    up, vp = get_velocity(V, t, setup)
    e = up .^ 2 .+ vp .^ 2
    reshape(abs.(fft(e)[k.+1, k.+1]), :)
end
espec = Figure()
ax =
    Axis(espec[1, 1]; xlabel = L"k", ylabel = L"\hat{e}(k)", xscale = log10, yscale = log10)
## ylims!(ax, (1e-20, 1))
scatter!(ax, kk, ehat; label = "Kinetic energy")
krange = LinRange(extrema(kk)..., 100)
lines!(ax, krange, 1e7 * krange .^ (-3); label = L"k^{-3}", color = :red)
axislegend(ax)
espec

# Solve unsteady problem
V, p = solve(V₀, p₀, (0.0, 1.0); setup, Δt = 0.0005, processors);

# Real time plot
rtp

# Energy history
ehist

# Energy spectrum
espec

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

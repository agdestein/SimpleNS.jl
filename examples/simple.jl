# Little LSP hack to get function signatures, go    #src
# to definition etc.                                #src
if isdefined(@__MODULE__, :LanguageServer)          #src
    include("../src/IncompressibleNavierStokes.jl") #src
    using .IncompressibleNavierStokes               #src
end                                                 #src

using FFTW
using GLMakie
using LinearAlgebra
using SparseArrays

function f(U, p)
    (ν, L) = p
    u, v, p = eachslice(u; dims = 3)
    nx, ny = size(u)
    Δx = L / nx
    Δy = L / ny

    up0 = circshift(u, (-1, 0))
    um0 = circshift(u, (1, 0))
    u0p = circshift(u, (0, 1))
    u0m = circshift(u, (0, -1))

    vp0 = circshift(v, (-1, 0))
    vm0 = circshift(v, (1, 0))
    v0p = circshift(v, (0, 1))
    v0m = circshift(v, (0, -1))

    du = @. ((u + up0)^2/4 - (um0 + u)^2 / 4) / Δx + ((u + u0p) / 2 * (v))
end

# A 2D grid is a Cartesian product of two vectors
n = 200
x = LinRange(0.0, 1.0, n + 1)
y = LinRange(0.0, 1.0, n + 1)

# Build setup and assemble operators
setup = Setup(x, y; viscosity_model);

# Initial conditions
K = n ÷ 2
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
Δp = pressure_poisson(pressure_solver, f);
V .-= (G * Δp ./  Ω);
p = pressure_additional_solve(pressure_solver, V, p, 0.0, setup);

V₀, p₀ = V, p

# Time interval
t_start, t_end = tlims = (0.0, 0.2)

# Iteration processors
logger = Logger()
observer = StateObserver(1, V₀, p₀, t_start)
writer = VTKWriter(; nupdate = 10, dir = "output/$name", filename = "solution")
## processors = [logger, observer, writer]
# processors = [logger, observer]
processors = [logger]

# Real time plot
rtp = real_time_plot(observer, setup)

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
problem = UnsteadyProblem(setup, V₀, p₀, (0.0, 0.1));
# problem = UnsteadyProblem(setup, V, p, tlims);
@time solve(
    problem,
    RK44();
    Δt = 0.001,
    processors,
    # pressure_solver = DirectPressureSolver(setup), # 4.2
    # pressure_solver = CGPressureSolver(setup;
    #     abstol = 10^-10,
    #     reltol = 10^-8,
    #     maxiter = 10,
    # ), # 3.0
    pressure_solver = FourierPressureSolver(setup), # 2.44
);

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
save_vtk(V, p, t_end, setup, "output/solution")

# Plot pressure
plot_pressure(setup, p)

# Plot velocity
plot_velocity(setup, V, t_end)

# Plot vorticity
plot_vorticity(setup, V, t_end)

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

# Device
# device = identity # For CPU
device = cu # For GPU

# Floating point type
T = Float32

# Domain size
lims = (T(0), T(1))

# Viscosity
Re = T(1000)
viscosity = 1 / Re

# DNS setup
N = 2^10
setup = get_setup(N, lims; viscosity);

# LES setup
M = 2^6
setup_les = get_setup(M, lims; viscosity);

# Create filtered training data
train = create_data(T(0.01), (T(0.0), T(0.5)); setup = device(setup), M, Δt = T(5f-5), nsubstep = 10, device)

pbar = zeros(T, M^2)
o = StateObserver(1, train.Vbar[:, 1], pbar, T(0))
rtp = real_time_plot(o, setup_les; type = :velocity)

for i in axes(train.Vbar, 2)
    Vb = train.Vbar[:, i]
    Fb = train.Fbar[:, i]
    F = momentum(Vb, Vb, pbar, setup_les; nopressure = true)
    o.state[] = (F - Fb, pbar, T(0))
    sleep(0.005)
end

# Initial conditions
K = N ÷ 2
V₀, p₀ = random_field(
    setup,
    K;
    A = T(1e8),
    σ = T(30),
    ## σ = 10,
    s = 5,
)

# Iteration processors
logger = Logger()
observer = StateObserver(20, V₀, p₀, T(0))
# writer = VTKWriter(; nupdate = 50, dir = "output/$name", filename = "solution")
# processors = [logger, observer, writer]
# processors = [logger, writer]
processors = [logger, observer]
# processors = [logger]

# Real time plot
rtp = real_time_plot(observer, setup; type = heatmap)

# Solve unsteady problem
@time V, p = solve(V₀, p₀, (0.0f0, 0.1f0); setup, Δt = 0.0005f0, processors);

# GPU
@time V, p = solve(
    cu(V₀), cu(p₀),
    # V, p,
    (0.0f0, 1.0f0);
    setup = cu(setup),
    Δt = 0.00005f0,
    processors,
);

# Real time plot
rtp

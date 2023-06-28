# Little LSP hack to get function signatures, go    #src
# to definition etc.                                #src
if isdefined(@__MODULE__, :LanguageServer)          #src
    include("../src/IncompressibleNavierStokes.jl") #src
    using .IncompressibleNavierStokes               #src
end                                                 #src

using CUDA
using FFTW
#md using CairoMakie
using GLMakie #!md
using IncompressibleNavierStokes
using LaTeXStrings
using LinearAlgebra
using SparseArrays

# Floating point type
T = Float64

# A 2D grid is a Cartesian product of two vectors
N = 2^10
lims = (T(0), T(1))
M = 2^6

# Viscosity
Re = T(1000)
viscosity = 1 / Re

# Build setup and assemble operators
setup = get_setup(N, lims; viscosity);
setup2 = get_setup(M, lims; viscosity);

# Initial conditions
K = N ÷ 2
V₀, p₀ = random_field(
    setup,
    K;
    A = T(1e6),
    # σ = T(30),
    σ = T(10),
    ## σ = 10,
    s = 5,
)

Wp = T.(create_top_hat_pressure(N, M))
W = T.(create_top_hat_velocity(N, M))

# plotmat(W)

Vbar = W * V₀
pbar = Wp * p₀

norm(setup.operators.M * V₀)
norm(setup2.operators.M * Vbar)

# Iteration processors
observer = StateObserver(20, V₀, p₀, T(0))
observer2 = StateObserver(20, Vbar, pbar, T(0))

# Real time plot
rtp = real_time_plot(observer, setup; type = heatmap, fieldname = :vorticity)
rtp2 = real_time_plot(observer2, setup2; type = heatmap, fieldname = :vorticity)

mkdir("output")
save("output/fine.png", rtp)
save("output/coarse.png", rtp2)

rtp
rtp2

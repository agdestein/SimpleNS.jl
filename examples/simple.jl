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

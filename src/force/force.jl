abstract type AbstractBodyForce{T} end

"""
    SteadyBodyForce{F1,F2,T}

Steady (constant) body force.
"""
struct SteadyBodyForce{F1,F2,T} <: AbstractBodyForce{T}
    fu::F1
    fv::F2
    F::Vector{T} # For storing discrete body force
end

"""
    SteadyBodyForce(fu, fv, grid)

Steady body force `f(x, y) = [fu(x, y), fv(x, y)]`. 
"""
function SteadyBodyForce(fu, fv, grid::Grid{T}) where {T}
    (; NV, indu, indv, xu, yu, xv, yv) = grid
    F = zeros(T, NV)
    F[indu] .= reshape(fu.(xu, yu), :)
    F[indv] .= reshape(fv.(xv, yv), :)
    SteadyBodyForce(fu, fv, F)
end

"""
    UnsteadyBodyForce{F1,F2,T}

Unsteady (time-dependent) body force.
"""
Base.@kwdef mutable struct UnsteadyBodyForce{F1,F2,T} <: AbstractBodyForce{T}
    fu::F1
    fv::F2
    F::Vector{T} # For storing discrete body force
end

"""
    UnsteadyBodyForce(fu, fv, grid)

Two-dimensional unsteady body force `f(x, y, t) = [fu(x, y, t), fv(x, y, t)]`. 
"""
function UnsteadyBodyForce(fu, fv, grid::Grid{T}) where {T}
    (; NV) = grid
    F = zeros(T, NV)
    UnteadyBodyForce(fu, fv, F)
end

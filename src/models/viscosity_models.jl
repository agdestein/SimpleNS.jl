"""
    AbstractViscosityModel

Abstract viscosity model.
"""
abstract type AbstractViscosityModel{T} end

"""
    LaminarModel(Re)

Laminar model with Reynolds number `Re`.
"""
Base.@kwdef struct LaminarModel{T} <: AbstractViscosityModel{T}
    Re::T # Reynolds number
end

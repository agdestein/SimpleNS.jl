"""
    AbstractConvectionModel

Abstract convection model.
"""
abstract type AbstractConvectionModel end

"""
    NoRegConvection()

Unregularized convection model.
"""
struct NoRegConvectionModel <: AbstractConvectionModel end

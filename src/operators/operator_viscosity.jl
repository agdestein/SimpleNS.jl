"""
    operator_viscosity(model, grid, boundary_conditions)

Classical turbulence modelling via the diffusive term
"""
function operator_viscosity end

operator_viscosity(::LaminarModel, grid, boundary_conditions) = (;)

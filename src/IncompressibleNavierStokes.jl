"""
    IncompressibleNavierStokes

Energy-conserving solvers for the incompressible Navier-Stokes equations.
"""
module IncompressibleNavierStokes

using FFTW
using Interpolations
using IterativeSolvers
using LinearAlgebra
using Printf
using SparseArrays
using Statistics
using UnPack
using WriteVTK: CollectionFile, paraview_collection, vtk_grid, vtk_save
using Makie

# Convenience notation
const ⊗ = kron

# Grid
include("grid/grid.jl")
include("grid/stretched_grid.jl")
include("grid/cosine_grid.jl")
include("grid/max_size.jl")

# Models
include("models/viscosity_models.jl")

# Types
include("force/force.jl")
include("boundary_conditions/boundary_conditions.jl")
include("operators/operators.jl")
include("setup.jl")

# Boundary condtions
include("boundary_conditions/bc_diff_stag.jl")
include("boundary_conditions/bc_general.jl")
include("boundary_conditions/bc_general_stag.jl")
include("boundary_conditions/get_bc_vectors.jl")

# Operators
include("operators/operator_averaging.jl")
include("operators/operator_convection_diffusion.jl")
include("operators/operator_divergence.jl")
include("operators/operator_interpolation.jl")
include("operators/operator_postprocessing.jl")
include("operators/operator_regularization.jl")

# Pressure solvers
include("solvers/pressure/pressure_solvers.jl")
include("solvers/pressure/pressure_poisson.jl")
include("solvers/pressure/pressure_additional_solve.jl")

# Time steppers
include("time_steppers/methods.jl")
include("time_steppers/tableaux.jl")
include("time_steppers/nstage.jl")
include("time_steppers/time_steppers.jl")
include("time_steppers/change_time_stepper.jl")
include("time_steppers/step.jl")
include("time_steppers/isexplicit.jl")
include("time_steppers/needs_startup_method.jl")
include("time_steppers/lambda_max.jl")

# Preprocess
include("preprocess/create_initial_conditions.jl")

# Processors
include("processors/processors.jl")
include("processors/initialize.jl")
include("processors/process.jl")
include("processors/finalize.jl")
include("processors/real_time_plot.jl")

# Momentum equation
include("momentum/bodyforce.jl")
include("momentum/compute_conservation.jl")
include("momentum/convection_components.jl")
include("momentum/convection.jl")
include("momentum/diffusion.jl")
include("momentum/momentum.jl")

# Problems
include("problems/problems.jl")
include("problems/is_steady.jl")

# Solvers
include("solvers/get_timestep.jl")
include("solvers/solve.jl")

# Utils
include("utils/get_lims.jl")

# Postprocess
include("postprocess/get_velocity.jl")
include("postprocess/get_vorticity.jl")
include("postprocess/get_streamfunction.jl")
include("postprocess/plot_force.jl")
include("postprocess/plot_grid.jl")
include("postprocess/plot_pressure.jl")
include("postprocess/plot_velocity.jl")
include("postprocess/plot_vorticity.jl")
include("postprocess/plot_streamfunction.jl")
include("postprocess/save_vtk.jl")

# Reexport
export @pack!

# Force
export SteadyBodyForce, UnsteadyBodyForce

# Models
export LaminarModel

# Processors
export AbstractProcessor, Logger, StateObserver, VTKWriter
export initialize!, process!, finalize!
export real_time_plot

# Setup
export Grid, Operators, BoundaryConditions, Setup

# 1D grids
export stretched_grid, cosine_grid

# Pressure solvers
export DirectPressureSolver, CGPressureSolver, FourierPressureSolver
export pressure_poisson, pressure_additional_solve

# Problems
export SteadyStateProblem, UnsteadyProblem, is_steady
export solve
export momentum

export create_initial_conditions, get_bc_vectors, get_velocity

export plot_force,
    plot_grid,
    plot_pressure,
    plot_streamfunction,
    plot_velocity,
    plot_vorticity,
    save_vtk

# ODE methods

export AdamsBashforthCrankNicolsonMethod, OneLegMethod

# Runge Kutta methods
export ExplicitRungeKuttaMethod, ImplicitRungeKuttaMethod, runge_kutta_method

# Explicit Methods
export FE11, SSP22, SSP42, SSP33, SSP43, SSP104, rSSPs2, rSSPs3, Wray3, RK56, DOPRI6

# Implicit Methods
export BE11, SDIRK34, ISSPm2, ISSPs3

# Half explicit methods
export HEM3, HEM3BS, HEM5

# Classical Methods
export GL1, GL2, GL3, RIA1, RIA2, RIA3, RIIA1, RIIA2, RIIA3, LIIIA2, LIIIA3

# Chebyshev methods
export CHDIRK3, CHCONS3, CHC3, CHC5

# Miscellaneous Methods
export Mid22, MTE22, CN22, Heun33, RK33C2, RK33P2, RK44, RK44C2, RK44C23, RK44P2

# DSRK Methods
export DSso2, DSRK2, DSRK3

# "Non-SSP" Methods of Wong & Spiteri
export NSSP21, NSSP32, NSSP33, NSSP53

end

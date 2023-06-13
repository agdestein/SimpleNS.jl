```@meta
CurrentModule = IncompressibleNavierStokes
```

# API Reference

```@docs
IncompressibleNavierStokes
Setup
```

## Force

```@docs
SteadyBodyForce
UnsteadyBodyForce
```

## Grid

```@docs
Grid
cosine_grid
get_dimension
max_size
stretched_grid
```

## Visocosity Models

```@docs
AbstractViscosityModel
LaminarModel
```

## Momentum

```@docs
bodyforce
compute_conservation
convection
convection_components
diffusion
momentum
momentum_allstage
turbulent_K
turbulent_viscosity
```

## Operators

```@docs
Operators
operator_averaging
operator_convection_diffusion
operator_divergence
operator_interpolation
operator_postprocessing
operator_turbulent_diffusion
```

## Postprocess

```@docs
get_velocity
get_vorticity
vorticity!
plot_force
plot_grid
plot_pressure
plot_velocity
plot_vorticity
save_vtk
```

## Preprocess

```@docs
create_initial_conditions
```

## Problems

```@docs
SteadyStateProblem
UnsteadyProblem
is_steady
```

## Processors

```@docs
AbstractProcessor
Logger
VTKWriter
StateObserver
initialize!
process!
finalize!
real_time_plot
```

## Solvers

```@docs
get_timestep
solve
solve_animate
```

### Pressure solvers

```@docs
AbstractPressureSolver
DirectPressureSolver
CGPressureSolver
FourierPressureSolver
pressure_additional_solve
pressure_poisson
```

## Time steppers

```@docs
AbstractODEMethod
AbstractRungeKuttaMethod
AdamsBashforthCrankNicolsonMethod
OneLegMethod
ExplicitRungeKuttaMethod
ImplicitRungeKuttaMethod

TimeStepper

change_time_stepper
isexplicit
lambda_conv_max
lambda_diff_max
needs_startup_method
nstage
runge_kutta_method
step
```

## Utils

```@docs
get_lims
```

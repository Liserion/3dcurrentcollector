# Micro app: single spherical particle, driven by the macro app via the
# CentroidMultiApp in macro_top_1c.i.
#
# This is the configuration validated on the cluster 2026-07-13 (job picard5)
# with Picard fixed-point coupling. Three things here were wrong in the old
# micro.i and are REQUIRED for the coupling to work — do not regress:
#   1. coef (Ds) = 2.2e-4, not 1.0e-5. The old value gave a ~25000-time-unit
#      particle (25x too slow for 1C): the surface quenches, the reaction
#      collapses, and the macro ghost-reacts at full current.
#   2. K2 = 171.375, identical to the macro's GlobalParams. With the old
#      K2 = 2.286 the macro and micro solve DIFFERENT physics and the Picard
#      iteration has no fixed point to converge to.
#   3. Tight tolerances. nl_rel_tol = 1e-4 silently accepts wrong surface
#      concentrations that poison the transferred cs.
#
# The parent MultiApp must run with sub_cycling = false: on the cluster's old
# ~/MOOSE framework, sub_cycling = true wipes the sub-app state to zero on any
# re-execution (every Picard iteration, every timestep cutback) — that is the
# historical "dies all the time" failure. With sub_cycling = false the parent
# imposes its dt, so the TimeStepper below is only a fallback.

[Mesh]
  type = GeneratedMesh
  dim = 1
  xmax = 0.5
  nx = 25
  coord_type = RSPHERICAL
[]

[Problem]
  type = FEProblem
[]

[Variables]
  [./Cs]
    family = LAGRANGE
    order = FIRST
    initial_condition = 0.5
  [../]
[]

[Kernels]
  [./dcsdt]
    type = TimeDerivative
    variable = Cs
  [../]
  [./csdiff]
    type = CoefDiffusion
    variable = Cs
    coef = 2.2e-4  # corrected particle diffusivity (Safari&Delacourt nano pair); old 1.0e-5 is infeasible at 1C
  [../]
[]

[GlobalParams]
  Cm = 0.1714785651793526
  K2 = 171.375   # MUST match the macro app's K2 (macro_top_1c.i GlobalParams)
[]

[BCs]
  [./bv_flux]
    type = ParticleBVPostBCKernel
    variable = Cs
    boundary = right
    pps_c2 = c2_from_macro
    pps_phi1 = phi1_from_macro
    pps_phi2 = phi2_from_macro
    MateChoice = 4
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK

  petsc_options_iname = '-pc_type -ksp_gmres_restart -pc_factor_mat_solver_type'
  petsc_options_value = ' lu       1601               mumps'

  nl_rel_tol = 8.5e-08
  nl_abs_tol = 1.5e-07

  # Fallback only: with sub_cycling = false the parent app sets dt directly.
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 5.0e-5
    optimal_iterations = 5
    growth_factor = 1.1
    cutback_factor = 0.5
  [../]
  dtmax = 1.0
  end_time = 36000.0   # must be >= the parent's end_time
[]

[Outputs]
  execute_on = 'INITIAL TIMESTEP_END'
  print_linear_residuals = false
  console = false
  csv = false
  exodus = false
[]

[Postprocessors]
  [./socp]
    type = ElementAverageValue
    variable = Cs
    execute_on = 'TIMESTEP_END'
  [../]
  [./cs_surface]
    type = SideAverageValue
    variable = Cs
    boundary = right
  [../]

  [./c2_from_macro]
    type = Receiver
  [../]
  [./phi1_from_macro]
    type = Receiver
  [../]
  [./phi2_from_macro]
    type = Receiver
  [../]
[]

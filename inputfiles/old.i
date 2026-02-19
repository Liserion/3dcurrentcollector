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
  [../]
[]

[ICs]
  [./cs_ic]
    type = FunctionIC
    variable = Cs
    # Non-uniform initial profile (center higher than surface)
    function = '0.55 - 0.10*(x/0.5)'
  [../]
[]

[Kernels]
  [./dcsdt]
    type = TimeDerivative
    variable = Cs
  [../]
  [./diff]
    type = CoefDiffusion
    variable = Cs
    coef = 1.0e-5
  [../]
[]

[GlobalParams]
#  Kappa = 1.61358
  #Kappa = 0.0
#  Chi = 2.5
  Cm = 0.1714785651793526
#  Sigma = 4426.43
#  K = 174.728
  K2 = 2.286
  #######################################

############################
#Omega = 0.0837492
[]
[BCs]
  [./test_flux]
    type = NeumannBC
    variable = Cs
    boundary = right
    value = 1e-1
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
  solve_type = NEWTON
  #solve_type = NEWTON

  petsc_options_iname = '-pc_type -ksp_gmres_restart -pc_factor_mat_solver_type'
  petsc_options_value = ' lu       1601               mumps'

  nl_rel_tol = 8.5e-08
  nl_abs_tol = 1.5e-07

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 5.0e-5
    optimal_iterations = 5
    growth_factor = 1.1
    cutback_factor = 0.5
  [../]
  dtmax = 1.0
  end_time = 36000.0
  #automatic_scaling=true
[]

[Outputs]
  execute_on = 'INITIAL TIMESTEP_END'
  print_linear_residuals = false
  console = true
  csv = false
  exodus = false
  #interval = 5
[]

[Postprocessors]
  [./socp]
    type = ElementAverageValue
    variable = Cs
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./cs_surface]
    type = SideAverageValue
    variable = Cs
    boundary = right
    execute_on = 'INITIAL TIMESTEP_END'
  [../]

#  [./J_from_macro]
#    type = Receiver
#  [../]
  [./c2_from_macro]
    type = Receiver
  [../]
  [./phi1_from_macro]
    type = Receiver
  [../]
  [./phi2_from_macro]
    type = Receiver
  [../]

  # [./mem]
  #   type= MemoryUsage
  # [../]
  # [./dofs]
  #   type= NumDOFs
  # [../]
  # [./elements]
  #   type= NumElems
  # [../]
[]

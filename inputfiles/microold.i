[Mesh]
  type = GeneratedMesh
  dim = 1
  xmax = 0.5
  nx = 25
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
    type = MatDiffusion
    variable = Cs
    diffusivity = 1.0e-5
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
  #solve_type = NEWTON

  # petsc_options_iname = '-pc_type -ksp_gmres_restart -pc_factor_mat_solver_type'
  # petsc_options_value = ' lu       1601               mumps'

  # nl_rel_tol = 8.5e-08
  # nl_abs_tol = 1.5e-07

  # [./TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 5.0e-5
  #   optimal_iterations = 5
  #   growth_factor = 1.1
  #   cutback_factor = 0.5
  # [../]
  # dtmax = 1.0
  # end_time = 36000.0
  # #automatic_scaling=true
  
  # Time stepping
  start_time = 0.0
  end_time = 36000.0
  dt = 1.0

  # PETSC solver
  # petsc_options_iname = '-pc_type -pc_hypre_type'
  # petsc_options_value = 'hypre boomeramg'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist                 '
  # automatic_scaling = true
  line_search = 'none'

  # Convergence criterion/criteria
  nl_max_its = 20
  l_max_its = 30
  nl_rel_tol = 1.0e-7
  nl_abs_tol = 1.0e-7
[]

[Outputs]
  execute_on = 'INITIAL TIMESTEP_END'
  #print_linear_residuals = false
  #console = false
  #csv = false
  #exodus = false
  #interval = 5
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

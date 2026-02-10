[Mesh]
  type = FileMesh
  file = macro_in.e
  construct_side_list_from_node_list = true
[]

[GlobalParams]
  D = 75.0
  Cm = 0.1714785651793526
  Sigma = 44.2643
  K = 1.74728
  #K2 = 342.75
  K2 = 2.286
  MateChoice = 4 # For LiFePO4
  a = 6.0
  #I = 3.26073 # 2C-rate
  #Omega = 0.08
[]

#[MeshModifiers]
#  [./block1]
#    type = SubdomainBoundingBox
#    block_id = 1
#    bottom_left = ' 45.0 0.0 0.0'
#    top_right   = '105.0 0.0 0.0'
#  [../]
#[]

[Variables]
  [./ce]
    initial_condition = 0.0874891
  [../]
  [./phis]
    initial_condition = 160.523
  [../]
  [./phie]
    initial_condition = 1.0e-12
  [../]
[]
[AuxVariables]
  [./cs]
    family = LAGRANGE
    order = FIRST
    initial_condition = 0.5
    block = cathode
  [../]
  [./soc]
    family = LAGRANGE
    order = FIRST
    initial_condition = 0.5
    block = cathode
  [../]
#   [./J]
#     family = LAGRANGE
#     order = FIRST
#     initial_condition = 1.0e-6
#     block = cathode
#   [../]
  [./Damage]
    family = LAGRANGE
    order = FIRST
    initial_condition = 0.0
  [../]
  [./SigmaH]
    family = LAGRANGE
    order = FIRST
    initial_condition = 0.0
  [../]
  [./RealC]
    family = LAGRANGE
    order = FIRST
    initial_condition = 1.0e-12
  [../]
  [./RealPhi1]
    family = LAGRANGE
    order = FIRST
    initial_condition = 1.0e-12
  [../]
  [./RealPhi2]
    family = LAGRANGE
    order = FIRST
    initial_condition = 1.0e-12
  [../]
[]

[Kernels]
  [./dcdt_separator]
    type = TimeDerivative
    variable = ce
    block = block_0
  [../]
  [./cdiff_separator]
    type = SeparatorCeKernel
    variable = ce
    PhiE = phie
    eps = 1.0
    block = block_0
  [../]
  [./phi1_separator]
    type = SeparatorPhiSKernel
    variable = phis
    block = block_0
  [../]
  [./phi2_separator]
    type = SeparatorPhiEKernel
    variable = phie
    Ce =  ce
    eps = 1.0
    block = block_0
  [../]
  ###############################
  ### For cathode
  [./dcdt_cathode]
    type = CoefTimeDerivative
    variable = ce
    Coefficient = 0.2
    block = cathode
  [../]
  [./cdiff_cathode]
    type = CathodeCeKernel
    variable = ce
    PhiS = phis
    PhiE = phie
    Cs = cs
    eps = 0.2
    Damage = Damage
    SigmaH = SigmaH
    block = cathode
  [../]
  [./phi1_cathode]
    type = CathodePhiSKernel
    variable = phis
    Ce = ce
    PhiE = phie
    Cs = cs
    eps = 0.2
    Damage = Damage
    SigmaH = SigmaH
    block = cathode
  [../]
  [./phi2_cathode]
    type = CathodePhiEKernel
    variable = phie
    Ce = ce
    PhiS = phis
    Cs = cs
    eps = 0.2
    Damage = Damage
    SigmaH = SigmaH
    block = cathode
  [../]
[]

[AuxKernels]
  [./getC]
    type = GetRealValueAuxKernel
    variable = RealC
    dof = ce
    CoefFactor = 22860.0
  [../]
  [./getPhi1]
    type = GetRealValueAuxKernel
    variable = RealPhi1
    dof = phis
    CoefFactor = 0.025690705483640282
  [../]
  [./getPhi2]
    type = GetRealValueAuxKernel
    variable = RealPhi2
    dof = phie
    CoefFactor = -0.025690705483640282
  [../]
#   [./getJ]
#     type = GetCellLevelFluxAuxKernel
#     variable = J
#     Ce = ce
#     Cs = cs
#     PhiE = phie
#     PhiS = phis
#     SigmaH = SigmaH
#     block = cathode
#   [../]
[]

[BCs]
  [./flux_c]
    type = ConstFluxForCeBC
    variable = ce
    boundary = top
    I = 0.0372588
  [../]
  [./flux_phi1]
    type = ConstFluxForPhiSBC
    variable = phis
    boundary = cat_cc
    I = 0.0070659
  [../]
  [./flux_phi2]
    type = ConstFluxForPhiEBC
    variable = phie
    boundary = top
    I = 0.0372588
  [../]
#  [./PhiE]
#    type = DirichletBC
#    variable = phie
#    value = 0.0
#    boundary = cat_cc
#  [../]
  [./PhiS]
    type = DirichletBC
    variable = phis
    value = 0.0
    boundary = top
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

  automatic_scaling = true
  compute_scaling_once = false

  nl_rel_tol = 8.5e-08
  nl_abs_tol = 1.5e-07

  # picard_max_its = 80
  # picard_rel_tol = 6.5e-08
  # picard_abs_tol = 1.0e-07

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0
    optimal_iterations = 10
    growth_factor = 1.2
    cutback_factor = 0.5
  [../]
  dtmax = 1.0
  end_time = 36000.0
  # num_steps = 2

  steady_state_detection = true
  steady_state_start_time = 12.0
  steady_state_tolerance = 9e-09
[]

[Outputs]
  # [./my_checkpoint]
  #   type = Checkpoint
  #   num_files = 5
  #   interval = 5000
  # [../]
  csv = true
  exodus = true
  execute_on = 'TIMESTEP_END'
  #interval = 2
  console = true
  print_linear_residuals = true
[]

[Debug]
  show_var_residual_norms = true
[]

[Postprocessors]
  [./cellv]
    type = SideAverageValue
    variable = RealPhi2
    boundary = top
    execute_on = 'TIMESTEP_END'
  [../]
  [./soc]
    type = ElementAverageValue
    variable = soc
    block = cathode
  [../]
  [./soe]
    type = ElementAverageValue
    variable = ce
    block = cathode
  [../]
  [./soet]
    type = ElementAverageValue
    variable = ce
  [../]
#  [./damage]
#    type = ElementAverageValue
#    variable = Damage
#    block = cathode
#  [../]
#  [./hystress]
#    type = ElementAverageValue
#    variable = SigmaH
#    block = cathode
#  [../]
#  [./socmax]
#    type = ElementExtremeValue
#    variable = soc
#    value_type = max
#  [../]
#  [./socmin]
#    type = ElementExtremeValue
#    variable = soc
#    value_type = min
#  [../]
  [./dt]
    type = TimestepSize
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

###################################################
### For the two-level framework
###################################################
[MultiApps]
  [./micro]
    type = CentroidMultiApp
    app_type = babblerAPP
    use_displaced_mesh = false
    execute_on = 'TIMESTEP_END'
    # sub_cycling = true
    sub_cycling = true
    input_files =  micro.i
    block = cathode
  [../]
[]

 [Transfers]
  [./c2_to_micro]
    type = MultiAppVariableValueSamplePostprocessorTransfer
    multi_app = micro
    execute_on = SAME_AS_MULTIAPP
    direction = to_multiapp
    source_variable = ce
    postprocessor = c2_from_macro
  [../]
  [./phi1_to_micro]
    type = MultiAppVariableValueSamplePostprocessorTransfer
    multi_app = micro
    execute_on = SAME_AS_MULTIAPP
    direction = to_multiapp
    source_variable = phis
    postprocessor = phi1_from_macro
  [../]
  [./phi2_to_micro]
    type = MultiAppVariableValueSamplePostprocessorTransfer
    multi_app = micro
    execute_on = SAME_AS_MULTIAPP
    direction = to_multiapp
    source_variable = phie
    postprocessor = phi2_from_macro
  [../]
#  [./J_to_micro]
#      type = MultiAppVariableValueSamplePostprocessorTransfer
#      multi_app = micro
#      direction = to_multiapp
#      source_variable = J
#      postprocessor = J_from_macro
#      execute_on = SAME_AS_MULTIAPP
#    [../]
#   #####################################
#   ## Cs_surface to macro
  [./cs_from_micro]
    type = MultiAppPostprocessorInterpolationTransfer
    multi_app = micro
    execute_on = 'TIMESTEP_END'
    direction = from_multiapp
    variable = cs
    postprocessor = cs_surface
    displaced_source_mesh = false
    displaced_target_mesh = false
    num_points = 3
    power = 2
    radius = -1
  [../]
  ## Cs_surface to macro
  [./soc_from_micro]
    type = MultiAppPostprocessorInterpolationTransfer
    multi_app = micro
    execute_on = 'TIMESTEP_END'
    direction = from_multiapp
    variable = soc
    postprocessor = socp
    displaced_source_mesh = false
    displaced_target_mesh = false
    num_points = 3
    power = 2
    radius = -1
  [../]
#   ###########################
#   # [./sigmaH_from_micro]
#   #   type = MultiAppPostprocessorInterpolationTransfer
#   #   multi_app = micro
#   #   execute_on = 'TIMESTEP_END'
#   #   direction = from_multiapp
#   #   variable = SigmaH
#   #   postprocessor = SigmaH
#   #   displaced_source_mesh = false
#   #   displaced_target_mesh = false
#   #   num_points = 3
#   #   power = 2
#   #   radius = -1
#   # [../]
#   #########################
#  [./damage_from_micro]
#    type = MultiAppPostprocessorInterpolationTransfer
#    multi_app = micro
#    execute_on = 'TIMESTEP_END'
#    direction = from_multiapp
#    variable = Damage
#    postprocessor = Damage
#    displaced_source_mesh = false
#    displaced_target_mesh = false
#    num_points = 3
#    power = 2
#    radius = -1
#  [../]
[]

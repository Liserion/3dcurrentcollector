[Mesh]
  type = FileMesh
  file = macro_in.e
  construct_side_list_from_node_list = true
[]

[GlobalParams]
  D = 3500.0
  Cm = 0.1714785651793526
  Sigma = 4426.43
  K = 87.364
  #K2 = 342.75
  K2 = 171.375
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
    initial_condition = 133.4
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
    block = 'block_0'
  [../]
  [./cdiff_separator]
    type = SeparatorCeKernel
    variable = ce
    PhiE = phie
    eps = 1.0
    block = 'block_0'
  [../]
  [./phi1_separator]
    type = SeparatorPhiSKernel
    variable = phis
    # A real separator is an electronic insulator. Without this override the
    # kernel inherits the cathode's Sigma=4426 from GlobalParams, which turns
    # the separator solid phase + the phis=0 top anchor into an electron
    # bypass around the reaction (ARL report Table 1 requires i1.n = 0 here).
    Sigma = 0.01
    block = 'block_0'
  [../]
  [./phi2_separator]
    type = SeparatorPhiEKernel
    variable = phie
    Ce =  ce
    eps = 1.0
    block = 'block_0'
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
  # TRUE first-principles 1C (with the insulating separator above, all current
  # must intercalate): I_top = (1-eps)*V_cathode/(3600*A_top*(1-t0)).
  # MEASURED mesh geometry (VolumePostprocessor/AreaPostprocessor, 2026-07-09):
  #   A_top = 1.14982e5, A_catcc = 5.92040e5, V_cathode = 2.13510e7
  # (the previously commented catcc=2.705984e5 and vol=1.904049e7 were stale —
  #  from a different mesh; with the wrong 0.4249 ratio the solid BC injects
  #  2.19x too much current and the excess poisons the voltage once the
  #  separator is sealed. Balanced ratio for THIS mesh: A_top/A_catcc = 0.19421)
  # With V = 2.1351e7: I_top(1C) = 0.8*2.1351e7/(3600*1.14982e5*0.989) = 0.04172
  # NOTE: Armin's 0.0372588 was 1C for HIS mesh volume; for this mesh use 0.0417.
  # The old fast-cell I=1.108 was a compensation for the separator bypass
  # (effective 0.14C); with the bypass sealed it would mean a real ~27C.
  [./flux_c]
    type = ConstFluxForCeBC
    variable = ce
    boundary = top
    I = 0.04172     # true 1C for this mesh   (old bypass-era value: 1.108)
  [../]
  [./flux_phi1]
    type = ConstFluxForPhiSBC
    variable = phis
    boundary = cat_cc
    I = 0.0081025   # true 1C balanced: I_top * 0.19421   (old: 0.470)
  [../]
  [./flux_phi2]
    type = ConstFluxForPhiEBC
    variable = phie
    boundary = top
    I = 0.04172     # true 1C for this mesh   (old: 1.108)
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
  line_search = bt
  nl_max_its = 30

  petsc_options_iname = '-pc_type -ksp_gmres_restart -pc_factor_mat_solver_type'
  petsc_options_value = ' lu       1501                mumps'

  nl_rel_tol = 1e-04
  nl_abs_tol = 1e-05

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0e-4
    optimal_iterations = 5
    growth_factor = 1.5
    cutback_factor = 0.5
  [../]
  dtmax = 1.0
  end_time = 12000.0

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
  print_linear_residuals = false
  console = true
  #interval = 2
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
    app_type = babblerApp
    use_displaced_mesh = false
    execute_on = 'TIMESTEP_END'
    # sub_cycling = true
    sub_cycling = true
    input_files = micro_fast.i
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

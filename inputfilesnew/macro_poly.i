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
  K2 = 2.286
  MateChoice = 4 # For LiFePO4
  a = 6.0
[]

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
  # Volume-averaged solid concentration (replaces micro solve)
  [./cs_avg]
    initial_condition = 0.5
    block = cathode
  [../]
[]

[AuxVariables]
  # Surface concentration from polynomial approximation
  [./cs]
    family = LAGRANGE
    order = FIRST
    initial_condition = 0.5
    block = cathode
  [../]
  # soc is now just cs_avg itself
  [./soc]
    family = LAGRANGE
    order = FIRST
    initial_condition = 0.5
    block = cathode
  [../]
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
  # ---- Separator ----
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
    Ce = ce
    eps = 1.0
    block = block_0
  [../]

  # ---- Cathode (ce, phis, phie — same as original) ----
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

  # ---- Polynomial approximation: cs_avg evolution ----
  [./dcsavg_dt]
    type = TimeDerivative
    variable = cs_avg
    block = cathode
  [../]
  [./csavg_source]
    type = CathodeCSAvgKernel
    variable = cs_avg
    Ce = ce
    PhiS = phis
    PhiE = phie
    Cs = cs
    eps = 0.2
    Damage = Damage
    SigmaH = SigmaH
    block = cathode
  [../]
[]

[AuxKernels]
  # Polynomial approximation: cs_surface = cs_avg - j_n*Rs/(5*Ds)
  [./compute_cs_surface]
    type = PolynomialCsSurfaceAux
    variable = cs
    cs_avg = cs_avg
    Ce = ce
    PhiS = phis
    PhiE = phie
    Damage = Damage
    SigmaH = SigmaH
    eps = 0.2
    execute_on = 'TIMESTEP_BEGIN NONLINEAR'
    block = cathode
  [../]
  # soc = cs_avg (volume average)
  [./compute_soc]
    type = ParsedAux
    variable = soc
    expression = 'cs_avg'
    coupled_variables = 'cs_avg'
    block = cathode
    execute_on = 'TIMESTEP_END'
  [../]
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
[]

[BCs]
  [./flux_c]
    type = ConstFluxForCeBC
    variable = ce
    boundary = top
    I = 0.416996
  [../]
  [./flux_phi1]
    type = ConstFluxForPhiSBC
    variable = phis
    boundary = cat_cc
    I = 0.094383
  [../]
  [./flux_phi2]
    type = ConstFluxForPhiEBC
    variable = phie
    boundary = top
    I = 0.416996
  [../]
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
  solve_type = NEWTON
  line_search = bt

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type -snes_linesearch_type -snes_max_it -ksp_max_it'
  petsc_options_value = ' lu       mumps                      bt                    50           200'

  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-06
  nl_max_its = 50

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0e-6
    optimal_iterations = 8
    growth_factor = 1.1
    cutback_factor = 0.4
  [../]
  dtmin = 1.0e-12
  dtmax = 5.0
  end_time = 12000.0

  steady_state_detection = true
  steady_state_start_time = 12.0
  steady_state_tolerance = 9e-09
[]

[Outputs]
  csv = true
  exodus = true
  execute_on = 'TIMESTEP_END'
  print_linear_residuals = false
  console = true
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
  [./dt]
    type = TimestepSize
  [../]
[]

# No MultiApps or Transfers — polynomial approximation replaces the micro solve

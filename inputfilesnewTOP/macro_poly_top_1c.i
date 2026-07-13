# ============================================================================
# Subramanian polynomial-approximation variant of macro_top_1c.i
#
# Instead of solving the full 1D spherical diffusion PDE per macro element
# (CentroidMultiApp + micro_fast.i), the solid phase is reduced to the
# two-parameter parabolic-profile model of Subramanian et al. (2005):
#
#   d(cs_avg)/dt = -(3/Rs) * j_n              (CsAvgODEKernelPoly, a = 3/Rs = 6)
#   cs_surf      = cs_avg - j_n*Rs/(5*Ds)     (PolynomialCsClosureKernel)
#
# cs_avg AND cs (= surface concentration) are nonlinear Variables solved
# monolithically with ce/phis/phie — the closure is an algebraic equation in
# the same Newton system, so there is no sub-cycled MultiApp, no transfers,
# and no lagged AuxVariable coupling. The cathode transport kernels
# (*KernelPoly) use the substituted closure flux j_n = (cs_avg-cs)*5*Ds/Rs,
# so the Butler-Volmer exponential appears only in the closure equation.
# Same mesh, BCs, parameters and outputs as macro_top_1c.i for comparison.
#
# ----------------------------------------------------------------------------
# VALIDATED 1C CONFIGURATION (2026-07-11): with the corrected particle
# diffusivity Ds = 2.2e-4 (Safari & Delacourt 2011 nano-LFP pair), the
# sealed separator, measured-geometry balanced currents, and consistent ICs,
# this input delivers a complete true-1C discharge: 2.5 V cutoff at
# t = 1705 s, soc 0.5 -> 0.979, effective rate 1.011C (results/poly_true1c).
# Validity of the parabolic closure: |j_n| <= (1 - cs_avg) * 5*Ds/Rs
# = (1 - cs_avg) * 2.2e-3, comfortably above the 1C demand of 4.7e-5 —
# rates up to several C are representable before the closure cap binds.
# ----------------------------------------------------------------------------
# ============================================================================

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
  K2 = 171.375
  MateChoice = 4 # For LiFePO4
  a = 6.0        # = 3/Rs with Rs = 0.5
[]

[Variables]
  [./ce]
    initial_condition = 0.0874891
  [../]
  # Equilibrium-consistent ICs (eta = 0 at t=0): phis matches its top
  # Dirichlet BC, phie = -U_ocv(cs=0.5) = -133.4938 in the same gauge.
  # This avoids the huge inconsistent-IC transient of the original input,
  # which a tightly-converged monolithic solve cannot step over.
  [./phis]
    initial_condition = 133.37   # = U_ocv(cs=0.5): eta = 0 at t=0 in this gauge
  [../]
  [./phie]
    initial_condition = 1.0e-12
  [../]
  [./cs_avg]
    initial_condition = 0.5
    block = cathode
  [../]
  [./cs]
    initial_condition = 0.5
    block = cathode
  [../]
[]

[AuxVariables]
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
  # The cathode transport kernels use the substituted closure flux
  #   j_n = (cs_avg - cs) * 5*Ds/Rs   (linear, bounded)
  # so the Butler-Volmer exponential appears ONLY in the cs closure equation
  # and all Jacobians entering the preconditioner are exact.
  [./cdiff_cathode]
    type = CathodeCeKernelPoly
    variable = ce
    PhiE = phie
    Cs = cs
    CsAvg = cs_avg
    eps = 0.2
    Ds = 2.2e-4  # corrected: Safari&Delacourt D1=1.18e-18 m2/s with nano radius 36.5nm (was 1.0e-5, source-less)
    Rs = 0.5
    block = cathode
  [../]
  [./phi1_cathode]
    type = CathodePhiSKernelPoly
    variable = phis
    Cs = cs
    CsAvg = cs_avg
    eps = 0.2
    Ds = 2.2e-4  # corrected: Safari&Delacourt D1=1.18e-18 m2/s with nano radius 36.5nm (was 1.0e-5, source-less)
    Rs = 0.5
    block = cathode
  [../]
  [./phi2_cathode]
    type = CathodePhiEKernelPoly
    variable = phie
    Ce = ce
    Cs = cs
    CsAvg = cs_avg
    eps = 0.2
    Ds = 2.2e-4  # corrected: Safari&Delacourt D1=1.18e-18 m2/s with nano radius 36.5nm (was 1.0e-5, source-less)
    Rs = 0.5
    block = cathode
  [../]
  ###############################
  ### Subramanian solid-phase ODE (replaces the micro MultiApp)
  [./dcs_avg_dt]
    type = TimeDerivative
    variable = cs_avg
    block = cathode
  [../]
  [./cs_avg_src]
    type = CsAvgODEKernelPoly
    variable = cs_avg
    Cs = cs
    Ds = 2.2e-4  # corrected: Safari&Delacourt D1=1.18e-18 m2/s with nano radius 36.5nm (was 1.0e-5, source-less)
    Rs = 0.5
    block = cathode
  [../]
  # Algebraic parabolic-profile closure: cs = cs_avg - j_n*Rs/(5*Ds),
  # solved implicitly in the same Newton system as all other variables.
  # This is the only equation containing the Butler-Volmer kinetics.
  [./cs_closure]
    type = PolynomialCsClosureKernel
    variable = cs
    cs_avg = cs_avg
    Ce = ce
    PhiS = phis
    PhiE = phie
    Ds = 2.2e-4  # corrected: Safari&Delacourt D1=1.18e-18 m2/s with nano radius 36.5nm (was 1.0e-5, source-less)
    Rs = 0.5
    Damage = Damage
    SigmaH = SigmaH
    block = cathode
  [../]
[]

[AuxKernels]
  # soc = particle-averaged concentration (same output field as MultiApp run)
  [./soc_from_cs_avg]
    type = GetRealValueAuxKernel
    variable = soc
    dof = cs_avg
    CoefFactor = 1.0
    block = cathode
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
  # MEASURED mesh geometry (2026-07-09): A_top = 1.14982e5,
  # A_catcc = 5.92040e5, V_cathode = 2.13510e7 -> balanced ratio
  # I_catcc = I_top * A_top/A_catcc = I_top * 0.19421.
  # First-principles 1C for this mesh: I_top = 0.04172.
  # Default here: I_top = 0.03 (~0.72C true) — validated, converges robustly,
  # intercalation = 100% of applied current with the insulating separator.
  # Higher currents exceed the parabolic closure cap (see header).
  [./flux_c]
    type = ConstFluxForCeBC
    variable = ce
    boundary = top
    I = 0.0417151   # true 1C (validated; crate_calc.py for other meshes/rates)
  [../]
  [./flux_phi1]
    type = ConstFluxForPhiSBC
    variable = phis
    boundary = cat_cc
    I = 0.0081016   # true 1C balanced: I_top * 0.19421 (measured areas)
  [../]
  [./flux_phi2]
    type = ConstFluxForPhiEBC
    variable = phie
    boundary = top
    I = 0.0417151   # true 1C (validated; crate_calc.py for other meshes/rates)
  [../]
  # Reference convention per the group's FE2 code (Betreuer, 2026-07-12):
  # electrolyte potential pinned to 0 at the current collector. phis then has
  # NO Dirichlet anywhere -> its equation is pure-Neumann and the injected
  # wire current MUST be consumed by the reaction (structurally leak-proof);
  # phis is anchored through the Butler-Volmer coupling to the pinned phie.
  [./PhiE_ref]
    type = DirichletBC
    variable = phie
    value = 0.0
    boundary = cat_cc
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

  # Tight tolerances are REQUIRED here: with the original loose values
  # (rel 1e-4 / abs 1e-5) the cs_avg/cs equations are accepted while still
  # unconverged and the soc trajectory is silently wrong by ~100x.
  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-06

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
  csv = true
  exodus = true
  execute_on = 'TIMESTEP_END'
  print_linear_residuals = false
  console = true
[]

[Postprocessors]
  # gauge-invariant cell voltage: V = RT/F*(phis@collector - phie@top)
  # (= RealPhi1@cat_cc + RealPhi2@top since RealPhi2 = -phie*RT/F)
  [./cellv_phis_cc]
    type = SideAverageValue
    variable = RealPhi1
    boundary = cat_cc
    execute_on = 'TIMESTEP_END'
  [../]
  [./cellv_phie_top]
    type = SideAverageValue
    variable = RealPhi2
    boundary = top
    execute_on = 'TIMESTEP_END'
  [../]
  [./cellv]
    type = ParsedPostprocessor
    function = 'cellv_phis_cc + cellv_phie_top'   # 'function' = old param name, works on all MOOSE versions
    pp_names = 'cellv_phis_cc cellv_phie_top'
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
  [./cs_surf_avg]
    type = ElementAverageValue
    variable = cs
    block = cathode
  [../]
  [./dt]
    type = TimestepSize
  [../]
[]

[Mesh]
  [./import]
    type = FileMeshGenerator
    file = 'wire_mesh.msh'
  []

  # Uniform scale: cathode z = 1.504 units → 170 µm, factor = 113.03
  [./scale]
    type = TransformGenerator
    input = import
    transform = SCALE
    vector_value = '113.03 113.03 113.03'
  []

  [./shift]
    type = TransformGenerator
    input = scale
    transform = TRANSLATE
    vector_value = '110.77 110.77 101.73'
  []

  [./catcc]
    type = SideSetsBetweenSubdomainsGenerator
    input = shift
    primary_block = 'cathode'
    paired_block = 'cat_inside'
    new_boundary = 'cat_cc'
  []

  [./strip_cat_inside]
    type = BlockDeletionGenerator
    input = catcc
    block = 'cat_inside'
  []

[]

[Variables]
  [./dummy]
  []
[]

[Kernels]
  [./time]
    type = TimeDerivative
    variable = dummy
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1e-6
  solve_type = PJFNK
[]

[Outputs]
  file_base = 'macro_in'
  [./exodus_out]
    type = Exodus
    execute_on = 'FINAL'
  []
  [./console_pps]
    type = Console
    execute_postprocessors_on = 'INITIAL'
  []
[]

# Print the geometry quantities needed to balance the BCs of macro.i.
# After running exodus.i, the console will show:
#   top_area, catcc_area, cathode_vol
# Then: I_top_target × A_top must equal I_phis_target × A_catcc.
[Postprocessors]
  [./top_area]
    type = AreaPostprocessor
    boundary = 'top'
    execute_on = 'INITIAL'
  [../]
  [./catcc_area]
    type = AreaPostprocessor
    boundary = 'cat_cc'
    execute_on = 'INITIAL'
  [../]
  [./cathode_vol]
    type = VolumePostprocessor
    block = 'cathode'
    execute_on = 'INITIAL'
  [../]
[]
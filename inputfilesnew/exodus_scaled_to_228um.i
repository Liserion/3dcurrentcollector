[Mesh]
  [./import]
    type = FileMeshGenerator
    file = 'firsttry.msh'
    uniform_refine = 0
  []

  # Optional: rescale Z so total thickness matches the original 228 µm cell thickness
  # New mesh bbox z_max = 298.782, so scale_z = 228/298.782 = 0.763098178605137
  [./scale_to_old_thickness]
    type = TransformGenerator
    input = import
    transform = SCALE
    vector_value = '1 1 0.763098178605137'   # 228 / 298.782
  []

  # 1) Create a sideset on the interface between cathode and cat_inside
  [./catcc]
    type = SideSetsBetweenSubdomainsGenerator
    input = scale_to_old_thickness
    primary_block = 'cathode'      # or the numeric ID, e.g. '1'
    paired_block  = 'cat_inside'   # or its ID, e.g. '2'
    new_boundary  = cat_cc
  []

  # 2) Now delete the cat_inside volume; the inner surface stays as cat_cc
  [./strip_cat_inside]
    type = BlockDeletionGenerator
    input = catcc
    block = 'cat_inside'           # or its ID, e.g. '2'
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
  [./exodus_out]
    type = Exodus
    execute_on = 'FINAL'
  []
[]

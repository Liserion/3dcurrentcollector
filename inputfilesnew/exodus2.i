[Mesh]
  [./import]
    type = FileMeshGenerator
    file = 'firsttry.msh'
  []

  # Uniformly scale the entire mesh so the overall size keeps its aspect ratio
  # while matching the original cell thickness (old z_max = 228 µm).
  # New mesh bbox z_max = 298.782 µm  =>  s = 228/298.782 = 0.763098178605137
  [./scale_uniform_to_old]
    type = TransformGenerator
    input = import
    transform = SCALE
    vector_value = '0.763098178605137 0.763098178605137 0.763098178605137'
  []

  # Create a sideset on the interface between cathode and cat_inside
  [./catcc]
    type = SideSetsBetweenSubdomainsGenerator
    input = scale_uniform_to_old
    primary_block = 'cathode'
    paired_block  = 'cat_inside'
    new_boundary  = cat_cc
  []

  # Delete the cat_inside volume; the inner surface stays as cat_cc
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
  [./exodus_out]
    type = Exodus
    execute_on = 'FINAL'
  []
[]

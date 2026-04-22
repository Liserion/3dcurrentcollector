[Mesh]
  [./import]
    type = FileMeshGenerator
    file = 'closed.msh'
  []

  # Scale to match old mesh dimensions (336x336x228 µm)
  [./scale]
    type = TransformGenerator
    input = import
    transform = SCALE
    vector_value = '172.3077 172.3077 131.4879'
  []

  [./catcc]
    type = SideSetsBetweenSubdomainsGenerator
    input = scale
    primary_block = 'cathode'
    paired_block  = 'cat_inside'
    new_boundary  = cat_cc
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
  [./exodus_out]
    type = Exodus
    execute_on = 'FINAL'
  []
[]
[Mesh]
  [./import]
    type = FileMeshGenerator
    file = 'closed.msh'
  []

  # 1) Create a sideset on the interface between cathode and cat_inside
  [./catcc]
    type = SideSetsBetweenSubdomainsGenerator
    input = import
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
[Mesh]
  [./import]
    type = FileMeshGenerator
    file = 'closed.msh'
  []

  # 1) Create the top sideset on block_0 at z=0.834
  [./top]
    type = SideSetsAroundSubdomainGenerator
    input = import
    block = 'block_0'
    new_boundary = 'top'
    normal = '0 0 1'
  []

  # 2) Create a sideset on the interface between cathode and cat_inside
  [./catcc]
    type = SideSetsBetweenSubdomainsGenerator
    input = top
    primary_block = 'cathode'
    paired_block  = 'cat_inside'
    new_boundary  = cat_cc
  []

  # 3) Delete the cat_inside volume; the inner surface stays as cat_cc
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
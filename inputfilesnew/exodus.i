[Mesh]
  type = FileMesh
  file = 'closed.msh'
[]

[Variables]
  [./dummy]
  [../]
[]

[Kernels]
  [./time]
    type = TimeDerivative
    variable = dummy
  [../]
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
  [../]
[]
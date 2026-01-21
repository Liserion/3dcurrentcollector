#include "SeparatorPhiSKernel.h"

registerMooseObject("babblerApp", SeparatorPhiSKernel);

InputParameters SeparatorPhiSKernel::validParams()
{
    InputParameters params = Kernel::validParams();
    params.addRequiredParam<Real>("Sigma", "conductivity of solid phase");
    return params;
}

SeparatorPhiSKernel::SeparatorPhiSKernel(const InputParameters &parameters)
: Kernel(parameters),
  _Sigma(getParam<Real>("Sigma"))
{
}

Real SeparatorPhiSKernel::computeQpResidual()
{
    return _Sigma * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real SeparatorPhiSKernel::computeQpJacobian()
{
    return _Sigma * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

#include "SeparatorPhiEKernel.h"

registerMooseObject("babblerApp", SeparatorPhiEKernel);

InputParameters SeparatorPhiEKernel::validParams()
{
    InputParameters params = Kernel::validParams();

    params.addRequiredParam<Real>("K", "conductivity of liquid phase");
    params.addRequiredParam<Real>("eps", "porosity");

    params.addRequiredCoupledVar("Ce", "concentration of electrolyte");

    return params;
}

SeparatorPhiEKernel::SeparatorPhiEKernel(const InputParameters &parameters)
: Kernel(parameters),
  _couple_c(coupledValue("Ce")),
  _grad_couple_c(coupledGradient("Ce")),
  _couple_c_var(coupled("Ce")),
  _K(getParam<Real>("K")),
  _eps(getParam<Real>("eps"))

{
    Keff = _K * _eps * sqrt(_eps);
}

Real SeparatorPhiEKernel::computeQpResidual()
{
    Real t0 = 0.0107907 + _couple_c[_qp] * 1.48837e-4;
    Keff = _K * _eps * sqrt(_eps);
    return Keff * (_grad_u[_qp] - (1 - t0) * _grad_couple_c[_qp] / _couple_c[_qp]) * _grad_test[_i][_qp];
}

Real SeparatorPhiEKernel::computeQpJacobian()
{
    Keff = _K * _eps * sqrt(_eps);
    return Keff * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

Real SeparatorPhiEKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
    Real t0 = 0.0107907 + _couple_c[_qp] * 1.48837e-4;
    Real dt0 = 1.48837e-4;
    Keff = _K * _eps * sqrt(_eps);
    
    if (jvar == _couple_c_var)
    {
        return Keff * dt0 * (_grad_couple_c[_qp] / _couple_c[_qp]) * _phi[_j][_qp] * _grad_test[_i][_qp]
               + Keff * (1 - t0) * (_grad_couple_c[_qp] / (_couple_c[_qp] * _couple_c[_qp])) * _phi[_j][_qp] * _grad_test[_i][_qp]
               - Keff * ((1 - t0) / _couple_c[_qp]) * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
    }

    return 0.0;
}

// Created by Armin on 29.10.2020
#include "SeparatorCeKernel.h"
#include "Moose.h"
registerMooseObject("babblerApp", SeparatorCeKernel);

InputParameters SeparatorCeKernel::validParams()
{
    InputParameters params = Kernel::validParams();

    params.addRequiredParam<Real>("D", "diffusivity");
    params.addRequiredParam<Real>("eps", "porosity");
    params.addRequiredParam<Real>("K", "conductivity");

    params.addRequiredCoupledVar("PhiE", "potential of electrolyte");

    return params;
}

SeparatorCeKernel::SeparatorCeKernel(const InputParameters &parameters)
: Kernel(parameters),
  _grad_couple_phi2(coupledGradient("PhiE")),
  _couple_phi2_var(coupled("PhiE")),
  _D(getParam<Real>("D")),
  _K(getParam<Real>("K")),
  _eps(getParam<Real>("eps"))
{
    Deff = _D * _eps;
    Keff = _K * _eps * sqrt(_eps);
}

Real SeparatorCeKernel::computeQpResidual()
{
    Real t0 = 0.0107907 + _u[_qp] * 1.48837e-4;
    Real dt0 = 1.48837e-4;

    return Deff * _grad_u[_qp] * _grad_test[_i][_qp]
           - dt0 * Keff * (_grad_couple_phi2[_qp] - (1 - t0) * _grad_u[_qp] / _u[_qp]) * _grad_u[_qp] * _test[_i][_qp];
}

Real SeparatorCeKernel::computeQpJacobian()
{
    Real t0 = 0.0107907 + _u[_qp] * 1.48837e-4;
    Real dt0 = 1.48837e-4;

    return Deff * _grad_phi[_j][_qp] * _grad_test[_i][_qp]
           - dt0 * Keff * dt0 * (_grad_u[_qp] / _u[_qp]) * _grad_u[_qp] * _phi[_j][_qp] * _test[_i][_qp]
           - dt0 * Keff * (1 - t0) * (_grad_u[_qp] * _grad_u[_qp] / (_u[_qp] * _u[_qp])) * _phi[_j][_qp] * _test[_i][_qp]
           + dt0 * Keff * (1 - t0) * (_grad_u[_qp] / _u[_qp]) * _grad_phi[_j][_qp] * _test[_i][_qp]
           - dt0 * Keff * (_grad_couple_phi2[_qp] - (1 - t0) * _grad_u[_qp] / _u[_qp]) * _grad_phi[_j][_qp] * _test[_i][_qp];
}

Real SeparatorCeKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
    Real t0 = 0.0107907 + _u[_qp] * 1.48837e-4;
    Real dt0 = 1.48837e-4;

    if (jvar == _couple_phi2_var)
    {
        return -dt0 * Keff * _grad_phi[_j][_qp] * _grad_u[_qp] * _test[_i][_qp];
    }

    return 0.0;
}

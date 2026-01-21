//created by Armin 29.10.2020

#include "ConstFluxForCeBC.h"

registerMooseObject("babblerApp", ConstFluxForCeBC);

InputParameters 
ConstFluxForCeBC::validParams()
{
    InputParameters params=IntegratedBCBase::validParams();

    params.addRequiredParam<Real>("I","current");
    params.addParam<Real>("ChargeTime",0.0,"0->forever,>=0 for charge time");

    return params;
}

ConstFluxForCeBC::ConstFluxForCeBC(const InputParameters &parameters)
:IntegratedBC(parameters),
_I(getParam<Real>("I")),
_ChargeTime(getParam<Real>("ChargeTime"))
{}

Real ConstFluxForCeBC::computeQpResidual()
{
    Real t0=0.0107907+_u[_qp]*1.48837e-4;

    if(_ChargeTime<=0.0)
    {
        return -_I*(1-t0)*_test[_i][_qp];
    }
    else
    {
        if(_t<=_ChargeTime)
        {
            return -_I*(1-t0)*_test[_i][_qp];
        }
        return 0.0;
    }
}

Real ConstFluxForCeBC::computeQpJacobian()
{
    const Real dt0=1.48837e-4;

    if(_ChargeTime>0.0)
    {
        return _I*dt0*_phi[_j][_qp]*_test[_i][_qp];
    }
    else
    {
        if(_t<=_ChargeTime)
        {
            return _I*dt0*_phi[_j][_qp]*_test[_i][_qp];
        }
        return 0.0;
    }

}

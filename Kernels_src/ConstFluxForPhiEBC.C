//created by Armin 29.10.2020

#include "ConstFluxForPhiEBC.h"


registerMooseObject("babblerApp", ConstFluxForPhiEBC);

InputParameters
ConstFluxForPhiEBC::validParams()
{
    InputParameters params=IntegratedBCBase::validParams();

    params.addRequiredParam<Real>("I","current");
    params.addParam<Real>("ChargeTime",0.0,"0->forever,>=0 for charge time");

    return params;
}

ConstFluxForPhiEBC::ConstFluxForPhiEBC(const InputParameters &parameters)
:IntegratedBC(parameters),
_I(getParam<Real>("I")),
_ChargeTime(getParam<Real>("ChargeTime"))
{}

Real ConstFluxForPhiEBC::computeQpResidual()
{
    if(_ChargeTime<=0.0)
    {
        return -_I*_test[_i][_qp];
    }
    else
    {
        if(_t<=_ChargeTime)
        {
            return -_I*_test[_i][_qp];
        }
        return 0.0;
    }
}

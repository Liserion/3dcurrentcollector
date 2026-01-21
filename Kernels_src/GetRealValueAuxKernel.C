
#include "GetRealValueAuxKernel.h"
registerMooseObject("babblerApp",GetRealValueAuxKernel);



InputParameters
GetRealValueAuxKernel::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("dof", "dof to be factored");
  params.addRequiredParam<Real>("CoefFactor", "factor coefficient");
  return params;
}

GetRealValueAuxKernel::GetRealValueAuxKernel(const InputParameters &parameters)
  : AuxKernel(parameters),
    _input_dof(coupledValue("dof")),
    _FactorValue(getParam<Real>("CoefFactor"))
{
}

Real
GetRealValueAuxKernel::computeValue()
{
  return _FactorValue * _input_dof[_qp];
}


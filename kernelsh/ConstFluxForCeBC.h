// created by Armin 29.10.2020

#pragma once

#include "IntegratedBC.h"

class ConstFluxForCeBC : public IntegratedBC
{
public:
    ConstFluxForCeBC(const InputParameters &parameters);

    static InputParameters validParams();

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    const Real &_I;
    const Real &_ChargeTime;
};

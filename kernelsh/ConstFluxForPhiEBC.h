// created by Armin 29.10.2020

#pragma once

#include "IntegratedBC.h"

class ConstFluxForPhiEBC : public IntegratedBC
{
public:
    ConstFluxForPhiEBC(const InputParameters &parameters);
    static InputParameters validParams();

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override { return 0.0; };

    const Real &_I;
    const Real &_ChargeTime;
};

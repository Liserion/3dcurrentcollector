#pragma once

#include "Kernel.h"

class SeparatorCeKernel : public Kernel
{
public:
    SeparatorCeKernel(const InputParameters &parameters);
    static InputParameters validParams();


protected:
    virtual Real computeQpResidual() override ;
    virtual Real computeQpJacobian() override ;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override ;

    const VariableGradient &_grad_couple_phi2;
    unsigned int _couple_phi2_var;

    const Real &_D,&_K,&_eps;

private:
    Real Deff,Keff;
};



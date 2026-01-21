#pragma once

#include "Kernel.h"

class SeparatorPhiEKernel : public Kernel
{
public:
    SeparatorPhiEKernel(const InputParameters &parameters);
    static InputParameters validParams();

protected:
    virtual Real computeQpResidual() override ;
    virtual Real computeQpJacobian() override ;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override ;

    const VariableValue &_couple_c;
    const VariableGradient &_grad_couple_c;
    unsigned int _couple_c_var;
    const Real &_K,&_eps;

private:
    Real Keff;
};



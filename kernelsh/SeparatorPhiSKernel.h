#pragma once

#include "Kernel.h"

class SeparatorPhiSKernel : public Kernel
{
public:
    SeparatorPhiSKernel(const InputParameters &parameters);
    static InputParameters validParams();

protected:
    virtual Real computeQpResidual() override ;
    virtual Real computeQpJacobian() override ;
    //virtual Real computeQpOffDiagJacobian(unsigned int jvar) override ;

    const Real &_Sigma;
};




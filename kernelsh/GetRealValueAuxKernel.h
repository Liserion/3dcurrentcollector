// Adapted to new format by [Your Name], [Date]
#pragma once

#include "AuxKernel.h"

class GetRealValueAuxKernel : public AuxKernel
{
public:
    static InputParameters validParams();
    explicit GetRealValueAuxKernel(const InputParameters &parameters);

protected:
    Real computeValue() override;

    // Member variables
    const VariableValue &_input_dof;
    const Real &_FactorValue;

    // Additional private members or utility functions can be declared below
};


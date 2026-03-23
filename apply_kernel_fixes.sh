#!/bin/bash
# Apply kernel fixes to babbler on HPC
# Run from: /work/scratch/jk98dago/babbler

BABBLER_DIR="$HOME/projects/babbler"
cd "$BABBLER_DIR" || { echo "ERROR: Cannot cd to $BABBLER_DIR"; exit 1; }

echo "=== Applying kernel fixes ==="

###############################################
# 1. ConstFluxForCeBC.h - add _RampTime member
###############################################
sed -i 's|    const Real &_ChargeTime;|    const Real \&_ChargeTime;\n    const Real \&_RampTime;|' include/kernels/ConstFluxForCeBC.h
echo "  [done] ConstFluxForCeBC.h"

###############################################
# 2. ConstFluxForPhiSBC.h - add _RampTime member
###############################################
sed -i 's|    const Real &_ChargeTime;|    const Real \&_ChargeTime;\n    const Real \&_RampTime;|' include/kernels/ConstFluxForPhiSBC.h
echo "  [done] ConstFluxForPhiSBC.h"

###############################################
# 3. ConstFluxForPhiEBC.h - add _RampTime member
###############################################
sed -i 's|    const Real &_ChargeTime;|    const Real \&_ChargeTime;\n    const Real \&_RampTime;|' include/kernels/ConstFluxForPhiEBC.h
echo "  [done] ConstFluxForPhiEBC.h"

###############################################
# 4. ConstFluxForCeBC.C - add RampTime param + ramp logic
###############################################
cat > src/Kernels/ConstFluxForCeBC.C << 'CEOF'
//created by Armin 29.10.2020

#include "ConstFluxForCeBC.h"

registerMooseObject("babblerApp", ConstFluxForCeBC);

InputParameters
ConstFluxForCeBC::validParams()
{
    InputParameters params=IntegratedBCBase::validParams();

    params.addRequiredParam<Real>("I","current");
    params.addParam<Real>("ChargeTime",0.0,"0->forever,>=0 for charge time");
    params.addParam<Real>("RampTime",0.0,"Ramp time: I scales as min(t/RampTime,1). 0=no ramp");

    return params;
}

ConstFluxForCeBC::ConstFluxForCeBC(const InputParameters &parameters)
:IntegratedBC(parameters),
_I(getParam<Real>("I")),
_ChargeTime(getParam<Real>("ChargeTime")),
_RampTime(getParam<Real>("RampTime"))
{}

Real ConstFluxForCeBC::computeQpResidual()
{
    Real t0=0.0107907+_u[_qp]*1.48837e-4;
    Real ramp = (_RampTime > 0.0) ? std::min(_t / _RampTime, 1.0) : 1.0;

    if(_ChargeTime<=0.0)
    {
        return -ramp*_I*(1-t0)*_test[_i][_qp];
    }
    else
    {
        if(_t<=_ChargeTime)
        {
            return -ramp*_I*(1-t0)*_test[_i][_qp];
        }
        return 0.0;
    }
}

Real ConstFluxForCeBC::computeQpJacobian()
{
    const Real dt0=1.48837e-4;
    Real ramp = (_RampTime > 0.0) ? std::min(_t / _RampTime, 1.0) : 1.0;

    if(_ChargeTime<=0.0)
    {
        return ramp*_I*dt0*_phi[_j][_qp]*_test[_i][_qp];
    }
    else
    {
        if(_t<=_ChargeTime)
        {
            return ramp*_I*dt0*_phi[_j][_qp]*_test[_i][_qp];
        }
        return 0.0;
    }

}
CEOF
echo "  [done] ConstFluxForCeBC.C"

###############################################
# 5. ConstFluxForPhiSBC.C - add RampTime param + ramp logic
###############################################
cat > src/Kernels/ConstFluxForPhiSBC.C << 'CEOF'
//created by Armin 29.10.2020

#include "ConstFluxForPhiSBC.h"


registerMooseObject("babblerApp", ConstFluxForPhiSBC);

InputParameters
ConstFluxForPhiSBC::validParams()
{
    InputParameters params=IntegratedBCBase::validParams();

    params.addRequiredParam<Real>("I","current");
    params.addParam<Real>("ChargeTime",0.0,"0->forever,>=0 for charge time");
    params.addParam<Real>("RampTime",0.0,"Ramp time: I scales as min(t/RampTime,1). 0=no ramp");

    return params;
}

ConstFluxForPhiSBC::ConstFluxForPhiSBC(const InputParameters &parameters)
:IntegratedBC(parameters),
_I(getParam<Real>("I")),
_ChargeTime(getParam<Real>("ChargeTime")),
_RampTime(getParam<Real>("RampTime"))
{}

Real ConstFluxForPhiSBC::computeQpResidual()
{
    Real ramp = (_RampTime > 0.0) ? std::min(_t / _RampTime, 1.0) : 1.0;

    if(_ChargeTime<=0.0)
    {
        return ramp*_I*_test[_i][_qp];
    }
    else
    {
        if(_t<=_ChargeTime)
        {
            return ramp*_I*_test[_i][_qp];
        }
        return 0.0;
    }
}
CEOF
echo "  [done] ConstFluxForPhiSBC.C"

###############################################
# 6. ConstFluxForPhiEBC.C - add RampTime param + ramp logic
###############################################
cat > src/Kernels/ConstFluxForPhiEBC.C << 'CEOF'
//created by Armin 29.10.2020

#include "ConstFluxForPhiEBC.h"


registerMooseObject("babblerApp", ConstFluxForPhiEBC);

InputParameters
ConstFluxForPhiEBC::validParams()
{
    InputParameters params=IntegratedBCBase::validParams();

    params.addRequiredParam<Real>("I","current");
    params.addParam<Real>("ChargeTime",0.0,"0->forever,>=0 for charge time");
    params.addParam<Real>("RampTime",0.0,"Ramp time: I scales as min(t/RampTime,1). 0=no ramp");

    return params;
}

ConstFluxForPhiEBC::ConstFluxForPhiEBC(const InputParameters &parameters)
:IntegratedBC(parameters),
_I(getParam<Real>("I")),
_ChargeTime(getParam<Real>("ChargeTime")),
_RampTime(getParam<Real>("RampTime"))
{}

Real ConstFluxForPhiEBC::computeQpResidual()
{
    Real ramp = (_RampTime > 0.0) ? std::min(_t / _RampTime, 1.0) : 1.0;

    if(_ChargeTime<=0.0)
    {
        return -ramp*_I*_test[_i][_qp];
    }
    else
    {
        if(_t<=_ChargeTime)
        {
            return -ramp*_I*_test[_i][_qp];
        }
        return 0.0;
    }
}
CEOF
echo "  [done] ConstFluxForPhiEBC.C"

echo ""
echo "=== All kernel fixes applied ==="
echo "Now run: make -j4"

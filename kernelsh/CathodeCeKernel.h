
#pragma once

#include "Kernel.h"

class CathodeCeKernel : public Kernel
{
public:
    static InputParameters validParams();
    explicit CathodeCeKernel(const InputParameters &parameters);

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
    const int &_MateChoice;
    const int &_IsDebug;
    const VariableValue &_couple_cs,&_couple_phi1,&_couple_phi2;
    const VariableGradient &_grad_couple_phi2;
    unsigned int _couple_phi1_var,_couple_phi2_var;


    const Real &_D,&_Cm,&_eps,&_K,&_K2,&_a;
    const Real _T;// gas constant, temperature and faraday constant


    // for damage  dependent && stress coupled over potential
    const VariableValue &_couple_damage,&_couple_sigmaH;
    const Real &_Omega;
private:
    Real Deff, Keff;
    Real Jeff, dJdc, dJdphi1, dJdphi2;

    Real OpenCircuitV(const int &matechoice, Real x);
    void BV(const Real &c, const Real &phi1, const Real &phi2, const Real &cs,
            Real &JEFF, Real &DJDC, Real &DJDPHI1, Real &DJDPHI2);
};

#pragma once

#include "Kernel.h"

class CathodePhiSKernel : public Kernel
{
public:
    CathodePhiSKernel(const InputParameters &parameters);
    static InputParameters validParams();

protected:
    virtual Real computeQpResidual() override ;
    virtual Real computeQpJacobian() override ;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override ;

    const VariableValue &_couple_c,&_couple_phi2,&_couple_cs;
    unsigned int _couple_c_var,_couple_phi2_var;


    const Real &_Sigma,&_eps,&_K2,&_Cm,&_a;
    const int &_MateChoice;
    const Real &_T;// gas constant, temperature and faraday constant


    // for damage  dependent && stress coupled over potential
    const VariableValue &_couple_damage,&_couple_sigmaH;
    const Real &_Omega;


private:
    Real Jeff,dJdc,dJdphi1,dJdphi2;

    Real OpenCircuitV(const int &matechoice,Real x);
    void BV(const Real &c,const Real &phi1,const Real &phi2,const Real &cs,
            Real &JEFF,Real &DJDC,Real &DJDPHI1,Real &DJDPHI2);
};




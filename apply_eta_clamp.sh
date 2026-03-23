#!/bin/bash
# Apply eta clamping to BV functions in cathode kernels and micro kernel
# This prevents exp() overflow at separator-cathode interface elements
# Run on HPC from any directory

BABBLER_DIR="$HOME/projects/babbler"
cd "$BABBLER_DIR" || { echo "ERROR: Cannot cd to $BABBLER_DIR"; exit 1; }

echo "=== Applying eta clamping fixes ==="

###############################################
# 1. CathodeCeKernel.C
###############################################
cat > src/Kernels/CathodeCeKernel.C << 'CEOF'
#include "CathodeCeKernel.h"

registerMooseObject("babblerApp", CathodeCeKernel);

InputParameters CathodeCeKernel::validParams()
{
    InputParameters params = Kernel::validParams();

    params.addRequiredParam<Real>("D", "diffusivity");
    params.addRequiredParam<Real>("Cm", "Max concentration of electrolyte");
    params.addRequiredParam<Real>("eps", "porosity");
    params.addRequiredParam<Real>("K", "conductivity of electrolyte");
    params.addRequiredParam<Real>("K2", "reaction rate");
    params.addParam<Real>("a", 3.0, "Area/Volume");
    params.addRequiredParam<int>("MateChoice","1--->TiS,"
                                               "2--->Mn2O4,"
                                               "3--->TiS modify,"
                                               "4--->LiFePO4,"
                                               "5--->LiFePO4 from Safari,"
                                               "6--->V2O5");

    params.addParam<int>("IsDebug", 0, "Debug flag for printing Js");
    params.addParam<Real>("T", 298.15, "temperature");

    params.addRequiredCoupledVar("PhiS", "potential for solid phase");
    params.addRequiredCoupledVar("PhiE", "potential for electrolyte phase");
    params.addRequiredCoupledVar("Cs", "surface concentration of solid particle");
    params.addRequiredCoupledVar("Damage", "damage from rve");
    params.addRequiredCoupledVar("SigmaH", "Hydrostatic stress from RVE homogenization");
    params.addParam<Real>("Omega", 0.0, "Partial molar volume of RVE material");

    return params;
}

CathodeCeKernel::CathodeCeKernel(const InputParameters & parameters) :
    Kernel(parameters),
    _MateChoice(getParam<int>("MateChoice")),
    _IsDebug(getParam<int>("IsDebug")),
    _couple_cs(coupledValue("Cs")),
    _couple_phi1(coupledValue("PhiS")),
    _couple_phi2(coupledValue("PhiE")),
    _grad_couple_phi2(coupledGradient("PhiE")),
    _couple_phi1_var(coupled("PhiS")),
    _couple_phi2_var(coupled("PhiE")),
    _D(getParam<Real>("D")),
    _Cm(getParam<Real>("Cm")),
    _eps(getParam<Real>("eps")),
    _K(getParam<Real>("K")),
    _K2(getParam<Real>("K2")),
    _a(getParam<Real>("a")),
    _T(getParam<Real>("T")),
    _couple_damage(coupledValue("Damage")),
    _couple_sigmaH(coupledValue("SigmaH")),
    _Omega(getParam<Real>("Omega"))
{
    // Constructor implementation
}

Real CathodeCeKernel::OpenCircuitV(const int &matechoice,Real x)
{
    Real V=0.0;

    const Real R=8.3144598;
    const Real F=96485.3329;

    if(matechoice==1)
    {
        V=2.17+(R*_T/F)*(log(fabs((1-x)/x))-16.2*x+8.1);
    }
    else if(matechoice==2)
    {
        V=4.06279+0.0677504*tanh(12.8268-21.8502*x)
          -0.105734*(pow(1.00167-x,-0.379571)-1.575994)
          -0.045*exp(-71.69*pow(x,8))
          +0.01*exp(-200.0*(x-0.19));
    }
    else if(matechoice==3)
    {
        V=2.17+R*_T*(-0.000558*x+8.10)/F;
    }
    else if(matechoice==4)
    {
        V=3.114559
          +4.438792*atan(-71.7352*x+70.85337)
          -4.240252*atan(-68.5605*x+67.730082);
    }
    else if(matechoice==5)
    {
        V=3.4324
          -0.8428*exp(-80.2493*pow(1-x,1.3198))
          -(3.2474e-6)*exp(20.2645*pow(1-x,3.8003))
          +(3.2482e-6)*exp(20.2646*pow(1-x,3.7995));
    }
    else if(matechoice==6)
    {
        V=3.3059
          +0.092769*tanh(-14.362*x+6.6874)
          -0.034252*exp(100*(x-0.96))
          +0.00724*exp(80.0*(0.01-x));
    }

    return V*F/(R*_T);
}

void CathodeCeKernel::BV(const Real &c, const Real &phi1, const Real &phi2,const Real &cs,
                        Real &JEFF,Real &DJDC,Real &DJDPHI1,Real &DJDPHI2)
{
    Real a=_a*(1-_eps);

    Real eta_raw=phi1-phi2-OpenCircuitV(_MateChoice,cs)-_Omega*_couple_sigmaH[_qp];

    // Clamp eta to prevent exp() overflow at interface elements
    // where phis transitions from separator (0) to cathode (133.4)
    bool eta_clamped = (eta_raw > 20.0 || eta_raw < -20.0);
    eta = std::max(std::min(eta_raw, 20.0), -20.0);

    // Clamp c to prevent division by zero and sqrt(negative)
    Real c_safe = std::max(c, 1.0e-10);
    Real sq = sqrt(std::max(_Cm*c_safe - c_safe*c_safe, 1.0e-20));

    JEFF=a*_K2*sq*(cs*exp(0.5*eta)-(1-cs)*exp(-0.5*eta));
    DJDC=a*0.5*_K2*(_Cm-2*c_safe)/sq*(cs*exp(0.5*eta)-(1-cs)*exp(-0.5*eta));
    if(eta_clamped)
    {
        // When eta is clamped, J is insensitive to phi1/phi2 changes
        DJDPHI1=0.0;
        DJDPHI2=0.0;
    }
    else
    {
        DJDPHI1=a*0.5*_K2*sq*(cs*exp(0.5*eta)+(1-cs)*exp(-0.5*eta));
        DJDPHI2=-a*0.5*_K2*sq*(cs*exp(0.5*eta)+(1-cs)*exp(-0.5*eta));
    }

    // add damage influence
    JEFF=(1-_couple_damage[_qp])*JEFF;
    DJDC=(1-_couple_damage[_qp])*DJDC;
    DJDPHI1=(1-_couple_damage[_qp])*DJDPHI1;
    DJDPHI2=(1-_couple_damage[_qp])*DJDPHI2;
    if(_couple_damage[_qp]>=0.9)
    {
        JEFF=0.0;
        DJDC=0.0;
        DJDPHI1=0.0;
        DJDPHI2=0.0;
    }
}


Real CathodeCeKernel::computeQpResidual()
{
    Real t0=0.0107907+_u[_qp]*1.48837e-4;
    Real dt0=1.48837e-4;

    Deff=_D*_eps;
    Keff=_K*_eps*sqrt(_eps);

    BV(_u[_qp],_couple_phi1[_qp],_couple_phi2[_qp],_couple_cs[_qp],Jeff,dJdc,dJdphi1,dJdphi2);

    if(_IsDebug)
    {
        std::cout<<"Js="<<Jeff
                 <<",dJdc="<<dJdc
                 <<",dJdphi1="<<dJdphi1
                 <<",dJdphi2="<<dJdphi2<<std::endl;
        std::cout<<"c="<<_u[_qp]
                 <<",Phi1="<<_couple_phi1[_qp]
                 <<",Phi2="<<_couple_phi2[_qp]
                 <<",Cs="<<_couple_cs[_qp]<<std::endl<<std::endl;
    }

    Real c_safe = std::max(_u[_qp], 1.0e-10);
    return Deff*_grad_u[_qp]*_grad_test[_i][_qp]
           -(1-t0)*Jeff*_test[_i][_qp]
           -dt0*Keff*(_grad_couple_phi2[_qp]-(1-t0)*_grad_u[_qp]/c_safe)*_grad_u[_qp]*_test[_i][_qp];
}

Real CathodeCeKernel::computeQpJacobian()
{
    Real t0=0.0107907+_u[_qp]*1.48837e-4;
    Real dt0=1.48837e-4;

    BV(_u[_qp],_couple_phi1[_qp],_couple_phi2[_qp],_couple_cs[_qp],Jeff,dJdc,dJdphi1,dJdphi2);

    Deff=_D*_eps;
    Keff=_K*_eps*sqrt(_eps);

    Real c_safe = std::max(_u[_qp], 1.0e-10);
    return Deff*_grad_phi[_j][_qp]*_grad_test[_i][_qp]
           -dt0*Keff*dt0*(_grad_u[_qp]/c_safe)*_grad_u[_qp]*_phi[_j][_qp]*_test[_i][_qp]
           -dt0*Keff*(1-t0)*(_grad_u[_qp]*_grad_u[_qp]/(c_safe*c_safe))*_phi[_j][_qp]*_test[_i][_qp]
           +dt0*Keff*(1-t0)*(_grad_u[_qp]/c_safe)*_grad_phi[_j][_qp]*_test[_i][_qp]
           -dt0*Keff*(_grad_couple_phi2[_qp]-(1-t0)*_grad_u[_qp]/c_safe)*_grad_phi[_j][_qp]*_test[_i][_qp]
           +dt0*Jeff*_phi[_j][_qp]*_test[_i][_qp]
           -(1-t0)*dJdc*_phi[_j][_qp]*_test[_i][_qp];
}


Real CathodeCeKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
    Real t0=0.0107907+_u[_qp]*1.48837e-4;
    Real dt0=1.48837e-4;



    Deff=_D*_eps;
    Keff=_K*_eps*sqrt(_eps);

    if(jvar==_couple_phi1_var)
    {
        BV(_u[_qp],_couple_phi1[_qp],_couple_phi2[_qp],_couple_cs[_qp],Jeff,dJdc,dJdphi1,dJdphi2);
        return -(1-t0)*dJdphi1*_phi[_j][_qp]*_test[_i][_qp];
    }
    else if(jvar==_couple_phi2_var)
    {
        BV(_u[_qp],_couple_phi1[_qp],_couple_phi2[_qp],_couple_cs[_qp],Jeff,dJdc,dJdphi1,dJdphi2);
        return -dt0*Keff*_grad_phi[_j][_qp]*_grad_u[_qp]*_test[_i][_qp]
               -(1-t0)*dJdphi2*_phi[_j][_qp]*_test[_i][_qp];
    }

    return 0.0;
}
CEOF
echo "  [done] CathodeCeKernel.C"

###############################################
# 2. CathodePhiSKernel.C
###############################################
cat > src/Kernels/CathodePhiSKernel.C << 'CEOF'
// Created by Armin on 29.10.2020
#include "CathodePhiSKernel.h"

registerMooseObject("babblerApp", CathodePhiSKernel);

InputParameters CathodePhiSKernel::validParams()
{
    InputParameters params = Kernel::validParams();

    params.addRequiredParam<Real>("Cm", "Max concentration of electrolyte");
    params.addRequiredParam<Real>("Sigma", "conductivity of solid phase");
    params.addRequiredParam<Real>("eps", "porosity");
    params.addRequiredParam<Real>("K2", "reaction rate");
    params.addParam<Real>("a", 3.0, "Area/Volume");

    params.addRequiredParam<int>("MateChoice","1--->TiS,"
                                              "2--->Mn2O4,"
                                              "3--->TiS modify,"
                                              "4--->LiFePO4,"
                                              "5--->LiFePO4 from Safari,"
                                              "6--->V2O5");

    // Physical constants
    params.addParam<Real>("T", 298.15, "temperature");

    params.addRequiredCoupledVar("Ce", "concentration of electrolyte");
    params.addRequiredCoupledVar("PhiE", "potential of electrolyte");
    params.addRequiredCoupledVar("Cs", "surface concentration of solid particle");
    params.addRequiredCoupledVar("Damage", "damage from rve");
    params.addRequiredCoupledVar("SigmaH", "Hydrostatic stress from RVE homogenization");

    params.addParam<Real>("Omega", 0.0, "Partial molar volume of RVE material");

    return params;
}

CathodePhiSKernel::CathodePhiSKernel(const InputParameters &parameters)
: Kernel(parameters),
  _couple_c(coupledValue("Ce")),
  _couple_phi2(coupledValue("PhiE")),
  _couple_cs(coupledValue("Cs")),
  _couple_c_var(coupled("Ce")),
  _couple_phi2_var(coupled("PhiE")),
  _Sigma(getParam<Real>("Sigma")),
  _eps(getParam<Real>("eps")),
  _K2(getParam<Real>("K2")),
  _Cm(getParam<Real>("Cm")),
 _a(getParam<Real>("a")),
  _MateChoice(getParam<int>("MateChoice")),
  _T(getParam<Real>("T")),
  _couple_damage(coupledValue("Damage")),
  _couple_sigmaH(coupledValue("SigmaH")),
  _Omega(getParam<Real>("Omega"))
{
  // Other constructor code
}

Real CathodePhiSKernel::OpenCircuitV(const int &matechoice,Real x)
{
    Real V=0.0;
    const Real R=8.3144598;
    const Real F=96485.3329;

    if(matechoice==1)
    {
        V=2.17+(R*_T/F)*(log(fabs((1-x)/x))-16.2*x+8.1);
    }
    else if(matechoice==2)
    {
        V=4.06279+0.0677504*tanh(12.8268-21.8502*x)
          -0.105734*(pow(1.00167-x,-0.379571)-1.575994)
          -0.045*exp(-71.69*pow(x,8))
          +0.01*exp(-200.0*(x-0.19));
    }
    else if(matechoice==3)
    {
        V=2.17+R*_T*(-0.000558*x+8.10)/F;
    }
    else if(matechoice==4)
    {
        V=3.114559
          +4.438792*atan(-71.7352*x+70.85337)
          -4.240252*atan(-68.5605*x+67.730082);
    }
    else if(matechoice==5)
    {
        V=3.4324
          -0.8428*exp(-80.2493*pow(1-x,1.3198))
          -(3.2474e-6)*exp(20.2645*pow(1-x,3.8003))
          +(3.2482e-6)*exp(20.2646*pow(1-x,3.7995));
    }
    else if(matechoice==6)
    {
        V=3.3059
          +0.092769*tanh(-14.362*x+6.6874)
          -0.034252*exp(100*(x-0.96))
          +0.00724*exp(80.0*(0.01-x));
    }

    return V*F/(R*_T);
}

void CathodePhiSKernel::BV(const Real &c, const Real &phi1, const Real &phi2,const Real &cs,
                            Real &JEFF,Real &DJDC,Real &DJDPHI1,Real &DJDPHI2)
{
    Real a=_a*(1-_eps);

    Real eta_raw=phi1-phi2-OpenCircuitV(_MateChoice,cs)-_Omega*_couple_sigmaH[_qp];

    // Clamp eta to prevent exp() overflow at interface elements
    bool eta_clamped = (eta_raw > 20.0 || eta_raw < -20.0);
    eta = std::max(std::min(eta_raw, 20.0), -20.0);

    Real c_safe = std::max(c, 1.0e-10);
    Real sq = sqrt(std::max(_Cm*c_safe - c_safe*c_safe, 1.0e-20));

    JEFF=a*_K2*sq*(cs*exp(0.5*eta)-(1-cs)*exp(-0.5*eta));
    DJDC=a*0.5*_K2*(_Cm-2*c_safe)/sq*(cs*exp(0.5*eta)-(1-cs)*exp(-0.5*eta));
    if(eta_clamped)
    {
        DJDPHI1=0.0;
        DJDPHI2=0.0;
    }
    else
    {
        DJDPHI1=a*0.5*_K2*sq*(cs*exp(0.5*eta)+(1-cs)*exp(-0.5*eta));
        DJDPHI2=-a*0.5*_K2*sq*(cs*exp(0.5*eta)+(1-cs)*exp(-0.5*eta));
    }

    // add damage influence
    JEFF=(1-_couple_damage[_qp])*JEFF;
    DJDC=(1-_couple_damage[_qp])*DJDC;
    DJDPHI1=(1-_couple_damage[_qp])*DJDPHI1;
    DJDPHI2=(1-_couple_damage[_qp])*DJDPHI2;
    if(_couple_damage[_qp]>=0.9)
    {
        JEFF=0.0;
        DJDC=0.0;
        DJDPHI1=0.0;
        DJDPHI2=0.0;
    }

}



Real CathodePhiSKernel::computeQpResidual()
{
    BV(_couple_c[_qp],_u[_qp],_couple_phi2[_qp],_couple_cs[_qp],Jeff,dJdc,dJdphi1,dJdphi2);

    return _Sigma*_grad_u[_qp]*_grad_test[_i][_qp]
           +Jeff*_test[_i][_qp];
}

Real CathodePhiSKernel::computeQpJacobian()
{
    BV(_couple_c[_qp],_u[_qp],_couple_phi2[_qp],_couple_cs[_qp],Jeff,dJdc,dJdphi1,dJdphi2);

    return _Sigma*_grad_phi[_j][_qp]*_grad_test[_i][_qp]
           +dJdphi1*_phi[_j][_qp]*_test[_i][_qp];
}

Real CathodePhiSKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
    if(jvar==_couple_c_var)
    {
        BV(_couple_c[_qp],_u[_qp],_couple_phi2[_qp],_couple_cs[_qp],Jeff,dJdc,dJdphi1,dJdphi2);
        return dJdc*_phi[_j][_qp]*_test[_i][_qp];
    }
    else if(jvar==_couple_phi2_var)
    {
        BV(_couple_c[_qp],_u[_qp],_couple_phi2[_qp],_couple_cs[_qp],Jeff,dJdc,dJdphi1,dJdphi2);
        return dJdphi2*_phi[_j][_qp]*_test[_i][_qp];
    }

    return 0.0;
}
CEOF
echo "  [done] CathodePhiSKernel.C"

###############################################
# 3. CathodePhiEKernel.C
###############################################
cat > src/Kernels/CathodePhiEKernel.C << 'CEOF'
#include "CathodePhiEKernel.h"


registerMooseObject("babblerApp",CathodePhiEKernel);

InputParameters CathodePhiEKernel::validParams()
{
    InputParameters params = Kernel::validParams();

    params.addRequiredParam<Real>("K","conductivity");
    params.addRequiredParam<Real>("Cm","Max concentration of electrolyte");
    params.addRequiredParam<Real>("eps","porosity");
    params.addRequiredParam<Real>("K2","reaction rate");
    params.addParam<Real>("a",3.0,"Area/Volume");

    params.addRequiredParam<int>("MateChoice","1--->TiS,"
                                              "2--->Mn2O4,"
                                              "3--->TiS modify,"
                                              "4--->LiFePO4,"
                                              "5--->LiFePO4 from Safari,"
                                              "6--->V2O5");

    // physical constant
    params.addParam<Real>("T",298.15,"temperature");


    params.addRequiredCoupledVar("Ce","concentration for electrolyte");
    params.addRequiredCoupledVar("PhiS","potential for solid phase");
    params.addRequiredCoupledVar("Cs","surface concentration of solid particle");

    params.addRequiredCoupledVar("Damage","damage from rve");
    params.addRequiredCoupledVar("SigmaH","Hydrostatic stress from RVE homogenezation");
    params.addParam<Real>("Omega",0.0,"Partial molar volume of RVE material");

    return params;
}

CathodePhiEKernel::CathodePhiEKernel(const InputParameters &parameters)
  : Kernel(parameters),
_couple_c(coupledValue("Ce")),
_couple_phi1(coupledValue("PhiS")),
_couple_cs(coupledValue("Cs")),
_grad_couple_c(coupledGradient("Ce")),
_couple_c_var(coupled("Ce")),
_couple_phi1_var(coupled("PhiS")),
    _K(getParam<Real>("K")),
_eps(getParam<Real>("eps")),
_K2(getParam<Real>("K2")),
    _Cm(getParam<Real>("Cm")),
_a(getParam<Real>("a")),
_MateChoice(getParam<int>("MateChoice")),
_T(getParam<Real>("T")),
_couple_damage(coupledValue("Damage")),
_couple_sigmaH(coupledValue("SigmaH")),
_Omega(getParam<Real>("Omega"))
{
    Keff=_K*_eps*sqrt(_eps);
}


Real CathodePhiEKernel::OpenCircuitV(const int &matechoice,Real x)
{
    Real V=0.0;

    const Real R=8.3144598;
    const Real F=96485.3329;

    if(matechoice==1)
    {
        V=2.17+(R*_T/F)*(log(fabs((1-x)/x))-16.2*x+8.1);
    }
    else if(matechoice==2)
    {
        V=4.06279+0.0677504*tanh(12.8268-21.8502*x)
          -0.105734*(pow(1.00167-x,-0.379571)-1.575994)
          -0.045*exp(-71.69*pow(x,8))
          +0.01*exp(-200.0*(x-0.19));
    }
    else if(matechoice==3)
    {
        V=2.17+R*_T*(-0.000558*x+8.10)/F;
    }
    else if(matechoice==4)
    {
        V=3.114559
          +4.438792*atan(-71.7352*x+70.85337)
          -4.240252*atan(-68.5605*x+67.730082);
    }
    else if(matechoice==5)
    {
        V=3.4324
          -0.8428*exp(-80.2493*pow(1-x,1.3198))
          -(3.2474e-6)*exp(20.2645*pow(1-x,3.8003))
          +(3.2482e-6)*exp(20.2646*pow(1-x,3.7995));
    }
    else if(matechoice==6)
    {
        V=3.3059
          +0.092769*tanh(-14.362*x+6.6874)
          -0.034252*exp(100*(x-0.96))
          +0.00724*exp(80.0*(0.01-x));
    }

    return V*F/(R*_T);
}

void CathodePhiEKernel::BV(const Real &c, const Real &phi1, const Real &phi2,const Real &cs,
                           Real &JEFF,Real &DJDC,Real &DJDPHI1,Real &DJDPHI2)
{
    Real a=_a*(1-_eps);

    Real eta_raw=phi1-phi2-OpenCircuitV(_MateChoice,cs)-_Omega*_couple_sigmaH[_qp];

    // Clamp eta to prevent exp() overflow at interface elements
    bool eta_clamped = (eta_raw > 20.0 || eta_raw < -20.0);
    eta = std::max(std::min(eta_raw, 20.0), -20.0);

    Real c_safe = std::max(c, 1.0e-10);
    Real sq = sqrt(std::max(_Cm*c_safe - c_safe*c_safe, 1.0e-20));

    JEFF=a*_K2*sq*(cs*exp(0.5*eta)-(1-cs)*exp(-0.5*eta));
    DJDC=a*0.5*_K2*(_Cm-2*c_safe)/sq*(cs*exp(0.5*eta)-(1-cs)*exp(-0.5*eta));
    if(eta_clamped)
    {
        DJDPHI1=0.0;
        DJDPHI2=0.0;
    }
    else
    {
        DJDPHI1=a*0.5*_K2*sq*(cs*exp(0.5*eta)+(1-cs)*exp(-0.5*eta));
        DJDPHI2=-a*0.5*_K2*sq*(cs*exp(0.5*eta)+(1-cs)*exp(-0.5*eta));
    }


    // add damage influence
    JEFF=(1-_couple_damage[_qp])*JEFF;
    DJDC=(1-_couple_damage[_qp])*DJDC;
    DJDPHI1=(1-_couple_damage[_qp])*DJDPHI1;
    DJDPHI2=(1-_couple_damage[_qp])*DJDPHI2;
    if(_couple_damage[_qp]>=0.9)
    {
        JEFF=0.0;
        DJDC=0.0;
        DJDPHI1=0.0;
        DJDPHI2=0.0;
    }

}


Real CathodePhiEKernel::computeQpResidual()
{
    BV(_couple_c[_qp],_couple_phi1[_qp],_u[_qp],_couple_cs[_qp],Jeff,dJdc,dJdphi1,dJdphi2);
    Real t0=0.0107907+_couple_c[_qp]*1.48837e-4;

    Keff=_K*_eps*sqrt(_eps);

    Real c_safe = std::max(_couple_c[_qp], 1.0e-10);
    return Keff*(_grad_u[_qp]-(1-t0)*_grad_couple_c[_qp]/c_safe)*_grad_test[_i][_qp]
           -Jeff*_test[_i][_qp];
}

Real CathodePhiEKernel::computeQpJacobian()
{
    BV(_couple_c[_qp],_couple_phi1[_qp],_u[_qp],_couple_cs[_qp],Jeff,dJdc,dJdphi1,dJdphi2);

    Keff=_K*_eps*sqrt(_eps);

    return Keff*_grad_phi[_j][_qp]*_grad_test[_i][_qp]
           -dJdphi2*_phi[_j][_qp]*_test[_i][_qp];
}

Real CathodePhiEKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
    Real t0=0.0107907+_couple_c[_qp]*1.48837e-4;
    Real dt0=1.48837e-4;

    Keff=_K*_eps*sqrt(_eps);


    if(jvar==_couple_c_var)
    {
        BV(_couple_c[_qp],_couple_phi1[_qp],_u[_qp],_couple_cs[_qp],Jeff,dJdc,dJdphi1,dJdphi2);
        Real c_safe = std::max(_couple_c[_qp], 1.0e-10);
        return Keff*dt0*(_grad_couple_c[_qp]/c_safe)*_phi[_j][_qp]*_grad_test[_i][_qp]
               +Keff*(1-t0)*(_grad_couple_c[_qp]/(c_safe*c_safe))*_phi[_j][_qp]*_grad_test[_i][_qp]
               -Keff*((1-t0)/c_safe)*_grad_phi[_j][_qp]*_grad_test[_i][_qp]
               -dJdc*_phi[_j][_qp]*_test[_i][_qp];
    }
    else if(jvar==_couple_phi1_var)
    {
        BV(_couple_c[_qp],_couple_phi1[_qp],_u[_qp],_couple_cs[_qp],Jeff,dJdc,dJdphi1,dJdphi2);
        return -dJdphi1*_phi[_j][_qp]*_test[_i][_qp];
    }
    return 0.0;
}
CEOF
echo "  [done] CathodePhiEKernel.C"

###############################################
# 4. ParticleBVPostBCKernel.C - eta clamping in micro BV
###############################################
cat > src/Kernels/ParticleBVPostBCKernel.C << 'CEOF'
#include "ParticleBVPostBCKernel.h"

registerMooseObject("babblerApp", ParticleBVPostBCKernel);

InputParameters
ParticleBVPostBCKernel::validParams()
{
    InputParameters params = IntegratedBC::validParams();

    params.addRequiredParam<PostprocessorName>("pps_c2", "name of postprocessor for electrolyte concentration c2");
    params.addRequiredParam<PostprocessorName>("pps_phi1", "name of postprocessor for solid phase potential phi1");
    params.addRequiredParam<PostprocessorName>("pps_phi2", "name of postprocessor for electrolyte phase potential phi2");

    params.addParam<Real>("Cm", 1.0, "Maximum concentration of solid particle");
    params.addParam<Real>("K2", 1.0, "Reaction rate for Butler-Volmer");
    params.addParam<Real>("T", 298.15, "Temperature (K)");
    params.addRequiredParam<int>("MateChoice",
                                 "1->TiS2, 2->Mn2O4, 3->TiS2 new, "
                                 "4->LiFePO4, 5->LiFePO4 from Safari, 6->V2O5");

    return params;
}

ParticleBVPostBCKernel::ParticleBVPostBCKernel(const InputParameters & parameters)
    : IntegratedBC(parameters),
      _pps_c2(getParam<PostprocessorName>("pps_c2")),
      _pps_phi1(getParam<PostprocessorName>("pps_phi1")),
      _pps_phi2(getParam<PostprocessorName>("pps_phi2")),
      _K2(getParam<Real>("K2")),
      _Cm(getParam<Real>("Cm")),
      _T(getParam<Real>("T")),
      _MateChoice(getParam<int>("MateChoice"))
{
}

void
ParticleBVPostBCKernel::OpenCircuitV(const Real & x, Real & u, Real & dudx)
{
    const Real R = 8.3144598, F = 96485.3329;

    if (_MateChoice == 1)
    {
        u = 2.17 + (R * _T / F) * (log(fabs((1 - x) / x)) - 16.2 * x + 8.1);
        dudx = (-16.2 * R * _T / F) * (x * x - x - 0.0617284) / (x * (x - 1));
    }
    else if (_MateChoice == 2)
    {
        u = 4.06279
            + 0.0677504 * tanh(12.8268 - 21.8502 * x)
            - 0.105734 * (pow(1.00167 - x, -0.379571) - 1.575994)
            - 0.045 * exp(-71.69 * pow(x, 8))
            + 0.01 * exp(-200.0 * (x - 0.19));

        dudx = -2.0 * exp(-200. * (x - 0.19))
               - 0.0401336 / pow(1.00167 - x, 1.37957)
               + 25.8084 * exp(-71.69 * pow(x, 8)) * pow(x, 7)
               - 1.48036 * Sech(12.8268 - 21.8502 * x) * Sech(12.8268 - 21.8502 * x);
    }
    else if (_MateChoice == 3)
    {
        u = 2.17 + R * _T * (-0.000558 * x + 8.10) / F;
        dudx = R * _T * -0.000558 / F;
    }
    else if (_MateChoice == 4)
    {
        // For LiFePO4
        u = 3.114559
            + 4.438792 * atan(-71.7352 * x + 70.85337)
            - 4.240252 * atan(-68.5605 * x + 67.730082);
        dudx = 4.438792 * (-71.7352) / (1 + pow(-71.7352 * x + 70.85337, 2))
             + (-4.240252) * (-68.5605) / (1 + pow(-68.5605 * x + 67.730082, 2));
    }
    else if (_MateChoice == 5)
    {
        u = 3.4324
            - 0.8428 * exp(-80.2493 * pow(1 - x, 1.3198))
            - (3.2474e-6) * exp(20.2645 * pow(1 - x, 3.8003))
            + (3.2482e-6) * exp(20.2646 * pow(1 - x, 3.7995));

        dudx = -89.2635 * exp(-80.2493 * pow(1 - x, 1.3198)) * pow(1 - x, 0.3198)
             + 0.000250086 * exp(20.2645 * pow(1 - x, 3.8003)) * pow(1 - x, 2.8003)
             - 0.000250104 * exp(20.2646 * pow(1 - x, 3.7995)) * pow(1 - x, 2.7995);
    }
    else if (_MateChoice == 6)
    {
        u = 3.3059
            + 0.092769 * tanh(-14.362 * x + 6.6874)
            - 0.034252 * exp(100 * (x - 0.96))
            + 0.00724 * exp(80.0 * (0.01 - x));

        dudx = 1.33235 * Sech(6.6874 - 14.362 * x) * Sech(6.6874 - 14.362 * x)
             - 3.4252 * exp(100 * (x - 0.96))
             - 0.5792 * exp(80 * (0.01 - x));
    }

    // Normalize to dimensionless units (divide by RT/F)
    u    = u * F / (R * _T);
    dudx = dudx * F / (R * _T);
}

void
ParticleBVPostBCKernel::BV(const Real & c, const Real & phi1, const Real & phi2,
                            const Real & cs, Real & J, Real & dJdc)
{
    Real U, dUdc;
    OpenCircuitV(cs, U, dUdc);

    Real eta_raw = phi1 - phi2 - U;

    // Clamp eta to prevent exp() overflow
    bool eta_clamped = (eta_raw > 20.0 || eta_raw < -20.0);
    Real eta = std::max(std::min(eta_raw, 20.0), -20.0);

    // Clamp c to prevent sqrt(negative) when electrolyte depletes
    Real c_safe = std::max(c, 1.0e-10);
    Real sq = sqrt(std::max((_Cm - c_safe) * c_safe, 1.0e-20));

    J = _K2 * sq * (cs * exp(0.5 * eta) - (1 - cs) * exp(-0.5 * eta));

    if(eta_clamped)
    {
        // When eta is clamped, deta/dcs = 0, so only direct cs derivatives remain
        dJdc = _K2 * sq * (exp(0.5 * eta) + exp(-0.5 * eta));
    }
    else
    {
        dJdc = _K2 * sq *
               (exp(0.5 * eta) - 0.5 * cs * dUdc * exp(0.5 * eta)
                + exp(-0.5 * eta) - 0.5 * (1 - cs) * dUdc * exp(-0.5 * eta));
    }
}

Real
ParticleBVPostBCKernel::computeQpResidual()
{
    _c2_value   = getPostprocessorValueByName(_pps_c2);
    _phi1_value = getPostprocessorValueByName(_pps_phi1);
    _phi2_value = getPostprocessorValueByName(_pps_phi2);

    BV(_c2_value, _phi1_value, _phi2_value, _u[_qp], J, dJdc);

    return J * _test[_i][_qp];
}

Real
ParticleBVPostBCKernel::computeQpJacobian()
{
    _c2_value   = getPostprocessorValueByName(_pps_c2);
    _phi1_value = getPostprocessorValueByName(_pps_phi1);
    _phi2_value = getPostprocessorValueByName(_pps_phi2);

    BV(_c2_value, _phi1_value, _phi2_value, _u[_qp], J, dJdc);

    return dJdc * _phi[_j][_qp] * _test[_i][_qp];
}
CEOF
echo "  [done] ParticleBVPostBCKernel.C"

echo ""
echo "=== All eta clamping fixes applied ==="
echo "Now run: make -j4"

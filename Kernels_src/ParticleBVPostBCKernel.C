

#include "ParticleBVPostBCKernel.h"

registerMooseObject("babblerApp", ParticleBVPostBCKernel);


InputParameters
ParticleBVPostBCKernel::validParams()
{
  InputParameters params = IntegratedBCBase::validParams();

  params.addRequiredParam<std::string>("pps_c2", "name of pps for c2");
  params.addRequiredParam<std::string>("pps_phi1", "name of pps for phi1");
  params.addRequiredParam<std::string>("pps_phi2", "name of pps for phi2");

  params.addRequiredParam<Real>("Cm", "The maximum concentration of electrolyte phase");
  params.addRequiredParam<Real>("K2", "The reaction rate for Bulter-Volmer reaction");
  params.addRequiredParam<int>("MateChoice",
                               "1->TiS2,"
                               "2->Mn2O4,"
                               "3->TiS2 new,"
                               "4->LiFePO4,"
                               "5->LiFePO4 from safari,"
                               "6->V2O5");

  params.addParam<Real>("T", 298.15, "Temperature(default 298.15)");

  return params;
}

ParticleBVPostBCKernel::ParticleBVPostBCKernel(const InputParameters &parameters)
:IntegratedBC(parameters),
  _pps_c2(getParam<std::string>("pps_c2")),
  _pps_phi1(getParam<std::string>("pps_phi1")),
  _pps_phi2(getParam<std::string>("pps_phi2")),
  _c2_value(getPostprocessorValueByName(_pps_c2)),
  _phi1_value(getPostprocessorValueByName(_pps_phi1)),
  _phi2_value(getPostprocessorValueByName(_pps_phi2)),
  _K2(getParam<Real>("K2")),
  _Cm(getParam<Real>("Cm")),
  _T(getParam<Real>("T")),
  _MateChoice(getParam<int>("MateChoice"))
{}

void ParticleBVPostBCKernel::OpenCircuitV(const Real &x, Real &u, Real &dudx)
{
    const Real R=8.3144598,F=96485.3329;


    if(_MateChoice==1)
    {
        // For TiS2
        // For LiyTiS2  (0.0<x<1.0)
        u=2.17+(R*_T/F)*(log(fabs((1-x)/x))-16.2*x+8.1);
        dudx=(-16.2*R*_T/F)*(x*x-x-0.0617284)/(x*(x-1));
    }
    else if(_MateChoice==2)
    {
        u=4.06279
          +0.0677504*tanh(12.8268-21.8502*x)
          -0.105734*(pow(1.00167-x,-0.379571)-1.575994)
          -0.045*exp(-71.69*pow(x,8))
          +0.01*exp(-200.0*(x-0.19));

        dudx=-2.0*exp(-200.*(x-0.19))
             -0.0401336/pow(1.00167-x,1.37957)
             +25.8084*exp(-71.69*pow(x,8))*pow(x,7)
             -1.48036*Sech(12.8268-21.8502*x)*Sech(12.8268-21.8502*x);
    }
    else if(_MateChoice==3)
    {
        u=2.17+R*_T*(-0.000558*x+8.10)/F;
        dudx=R*_T*-0.000558/F;
    }
    else if(_MateChoice==4)
    {
        // For LiFePO4
        // taken from "Discharge Model for the Lithium Iron-Phosphate Electrode"
        // Venkat Srinivasan, DOI: 10.1149/1.1785012
        u=3.114559
          +4.438792*atan(-71.7352*x+70.85337)
          -4.240252*atan(-68.5605*x+67.730082);
        dudx=4.438792*(-71.7352)/(1+pow(-71.7352*x+70.85337,2))
            +(-4.240252)*(-68.5605)/(1+pow(-68.5605*x+67.730082,2));
    }
    else if(_MateChoice==5)
    {
        // For LiFePO4
        // taken from Modeling of a Commercial Graphite/LiFePO4 Cell
        // Safari   doi: 10.1149/1.3567007
        u=3.4324
          -0.8428*exp(-80.2493*pow(1-x,1.3198))
          -(3.2474e-6)*exp(20.2645*pow(1-x,3.8003))
          +(3.2482e-6)*exp(20.2646*pow(1-x,3.7995));

        dudx=-89.2635*exp(-80.2493*pow(1-x,1.3198))*pow(1-x,0.3198)
            +0.000250086*exp(20.2645*pow(1-x,3.8003))*pow(1-x,2.8003)
            -0.000250104*exp(20.2646*pow(1-x,3.7995))*pow(1-x,2.7995);
    }
    else if(_MateChoice==6)
    {
        // For LiFePO4
        // taken from Modeling of a Commercial Graphite/LiFePO4 Cell
        // Safari   doi: 10.1149/1.3567007
        u=3.3059
          +0.092769*tanh(-14.362*x+6.6874)
          -0.034252*exp(100*(x-0.96))
          +0.00724*exp(80.0*(0.01-x));

        dudx=1.33235*Sech(6.6874-14.362*x)*Sech(6.6874-14.362*x)
            -3.4252*exp(100*(x-0.96))
            -0.5792*exp(80*(0.01-x));
    }

    u=u*F/(R*_T);
    dudx=dudx*F/(R*_T);
}

void ParticleBVPostBCKernel::BV(const Real &c, const Real &phi1, const Real &phi2, const Real &cs, Real &J, Real &dJdc)
{
    Real U,dUdc;

    OpenCircuitV(cs,U,dUdc);

    Real eta=phi1-phi2-U;

    J=_K2*sqrt((_Cm-c)*c)*(cs*exp(0.5*eta)-(1-cs)*exp(-0.5*eta));

    dJdc=_K2*sqrt((_Cm-c)*c)*
         (exp(0.5*eta)-0.5*cs*dUdc*exp(0.5*eta)
          +exp(-0.5*eta)-0.5*(1-cs)*dUdc*exp(-0.5*eta));
}


//*************************
Real ParticleBVPostBCKernel::computeQpResidual()
{
    Real c2 = _c2_value;
    Real phi1 = _phi1_value;
    Real phi2 = _phi2_value;

    BV(c2,phi1,phi2,_u[_qp],J,dJdc);

    return J*_test[_i][_qp];
}

Real ParticleBVPostBCKernel::computeQpJacobian()
{
    Real c2 = _c2_value;
    Real phi1 = _phi1_value;
    Real phi2 = _phi2_value;

    BV(c2,phi1,phi2,_u[_qp],J,dJdc);

    return dJdc*_phi[_j][_qp]*_test[_i][_qp];
}

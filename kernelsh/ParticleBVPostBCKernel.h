//created by Armin 29.10.2020

#pragma once




#include "IntegratedBC.h"

class ParticleBVPostBCKernel;

//InputParameters
//ParticleBVPostBCKernel::validParams()

class ParticleBVPostBCKernel:public IntegratedBC
{
public:
  ParticleBVPostBCKernel(const InputParameters & parameters);
  static InputParameters validParams();

protected:
    virtual Real computeQpResidual() override ;
    virtual Real computeQpJacobian() override ;
  PostprocessorName _pps_c2;
  PostprocessorName _pps_phi1;
  PostprocessorName _pps_phi2;
  Real _c2_value;
  Real _phi1_value;
  Real _phi2_value;
  Real _K2;
  Real _Cm;
  Real _T;
  int _MateChoice;
     Real J,dJdc;


     Real Sech(const Real &x)
    {
        return 2.0/(exp(x)+exp(-x));
    }

    void OpenCircuitV(const Real &x,Real &u,Real &dudx);
    void BV(const Real &c,const Real &phi1,const Real &phi2,
             const Real &cs,Real &J,Real &dJdc);
};

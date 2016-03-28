////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Feb 10 11:14:09 EST 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <cmath>

#include "TBField1DZRegularized.h"

int RK4Test ();
double By (double const Z, double const LambdaUndulator, double const PeakBy);
void rk4(float y[], float dydx[], int n, float x, float h, float yout[], void (*derivs)(float, float [], float []));
void derivs(float t, float y[], float dydt[]);


TBField1DZRegularized* TBF;

double const kPI = 3.14159265359;
double const k2PI = 2 * 3.14159265359;
double const BetaZ = 1.0 - 1E-16;
double const Beta = BetaZ;
double const Gamma = 1.0 / sqrt(1.0 - Beta*Beta);
double const kEMass = 9.10938356E-31;
double const kECharge = -1.60217662E-19;

double const kC_SI = 299792458.;
//double const XStart = -1.8076921081543;
double const XStart = -0.9;
double const XStop  =  0.9;


void derivs(float t, float x[], float dxdt[])
{
  dxdt[0] = x[1];
  dxdt[1] = -(kECharge * BetaZ * kC_SI) / (Gamma * kEMass) * TBF->GetByAtZ(XStart + t * BetaZ * kC_SI);
  std::cout << "By: " << XStart + t * BetaZ * kC_SI << "  " << TBF->GetByAtZ(XStart + t * BetaZ * kC_SI) << std::endl;

 
  return;
}



void rk4(float y[], float dydx[], int n, float x, float h, float yout[], void (*derivs)(float, float [], float []))
{
  int i;
  float xh,hh,h6,*dym,*dyt,*yt;
  dym=new float[n];
  dyt=new float[n];
  yt=new float[n];

  hh=h*0.5;
  h6=h/6.0;
  xh=x+hh;

  for (i=0;i<n;i++) {
    yt[i]=y[i]+hh*dydx[i];
  }
  (*derivs)(xh,yt,dyt);
  for (i=0;i<n;i++) {
    yt[i]=y[i]+hh*dyt[i];
  }
  (*derivs)(xh,yt,dym);
  for (i=0;i<n;i++) {
    yt[i]=y[i]+h*dym[i];
    dym[i] += dyt[i];
  }
  (*derivs)(x+h,yt,dyt);
  for (i=0;i<n;i++) {
    yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
  }
  
  delete [] dym;
  delete [] dyt;
  delete [] yt;

  return;
}


















int RK4Test ()
{


  int const N = 2;
  float x[N] = {0.0, 0.0};
  float dxdt[N];


  derivs(0, x, dxdt);

  int const NPoints = 10000;
  float const h = (XStop - XStart) / (BetaZ * kC_SI) / (NPoints - 1);
  std::cout << h << std::endl;

  for (int i = 0; i != NPoints; ++i) {
    float t = h * i;
    derivs(t, x, dxdt);
    rk4(x, dxdt, N, t, h, x, derivs);
    std::cout << "Position: " << XStart + kC_SI * BetaZ * (t+h) << "  " << x[0] << std::endl;
    std::cout << "dydx: " << i << " " << dxdt[1] << std::endl;
  }



  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InBFieldFile]" << std::endl;
    return 1;
  }


  TBF = new TBField1DZRegularized(argv[1]);
  TBF->SaveAs("del2.dat");
  std::cout << TBF->GetByAtZ(0.01) << std::endl;
  std::cout << TBF->GetByAtZ(0.11) << std::endl;
  std::cout << TBF->GetByAtZ(0.21) << std::endl;
  std::cout << TBF->GetByAtZ(0.31) << std::endl;
  std::cout << TBF->GetByAtZ(0.41) << std::endl;
  std::cout << TBF->GetByAtZ(0.49) << std::endl;
  std::cout << TBF->GetByAtZ(0.50) << std::endl;
  std::cout << TBF->GetByAtZ(0.40) << std::endl;


  RK4Test();

  delete TBF;

  return 0;
}

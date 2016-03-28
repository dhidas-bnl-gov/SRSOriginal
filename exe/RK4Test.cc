////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Feb 10 11:14:09 EST 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <cmath>

#include "TBField3DZRegularized.h"
#include "TVector3D.h"

#include "TGraph.h"
#include "TGraph2D.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFrame.h"

int RK4Test ();
double By (double const Z, double const LambdaUndulator, double const PeakBy);
void rk4(double y[], double dydx[], int n, double x, double h, double yout[], void (*derivs)(double, double [], double []));
void derivs(double t, double y[], double dydt[]);


TBField3DZRegularized* TBF;



double const EEnergy = 3.0;
double const EMass   = 0.000511;
double const EGamma  = EEnergy / EMass;
double const EBeta   = sqrt(1.0 - 1.0/(EGamma*EGamma));

double const kPI = 3.14159265359;
double const k2PI = 2 * 3.14159265359;
double const BetaZ = EBeta;
double const Beta = BetaZ;
double const Gamma = 1.0 / sqrt(1.0 - Beta*Beta);
double const kEMass = 9.10938356E-31;
double const kECharge = -1.60217662E-19;

double const BetaX = 0;
double const BetaY = 0;



double const kC_SI = 299792458.;
//double const XStart = -1.8076921081543;
double const XStart = -1.5;
double const XStop  =  1.5;








void derivs(double t, double x[], double dxdt[])
{
  //double b = sqrt(x[1]*x[1] + x[3]*x[3] + x[5]*x[5])/kC_SI;
  //double g = 1./sqrt(1-b);
  dxdt[0] = x[1];
  dxdt[1] = -(kECharge * x[5]) / (Gamma * kEMass) * TBF->GetBy(x[0], x[2], x[4]) + (kECharge * x[3]) / (Gamma * kEMass) * TBF->GetBz(x[0], x[2], x[4]);
  dxdt[2] = x[3];                                                                                                                            
  dxdt[3] =  (kECharge * x[5]) / (Gamma * kEMass) * TBF->GetBx(x[0], x[2], x[4]) - (kECharge * x[1]) / (Gamma * kEMass) * TBF->GetBz(x[0], x[2], x[4]);
  dxdt[4] = x[5];                                                                                                                            
  dxdt[5] =  (kECharge * x[1]) / (Gamma * kEMass) * TBF->GetBy(x[0], x[2], x[4]) - (kECharge * x[3]) / (Gamma * kEMass) * TBF->GetBx(x[0], x[2], x[4]);

  return;
}



void rk4(double y[], double dydx[], int n, double x, double h, double yout[], void (*derivs)(double, double [], double []))
{
  int i;
  double xh,hh,h6,*dym,*dyt,*yt;
  dym=new double[n];
  dyt=new double[n];
  yt=new double[n];

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
















double By (double const Z, double const LambdaUndulator, double const PeakBy)
{
  return PeakBy * cos(Z * k2PI / LambdaUndulator);
}



int RK4Test ()
{


  int const N = 6;
  double x[N] = {0.0, 0.0, 0.0, 0.0, (double) XStart, (double) (BetaZ * kC_SI)};
  double dxdt[N];


  int const NPoints = 10000;
  double const h = (XStop - XStart) / (BetaZ * kC_SI) / (NPoints - 1);

  // Graphs for viewing trajectory
  TGraph gXZ(NPoints);
  TGraph gYZ(NPoints);
  TGraph gYX(NPoints);

  TGraph gVX(NPoints);
  TGraph gVY(NPoints);
  TGraph gVZ(NPoints);

  TGraph gAX(NPoints);
  TGraph gAY(NPoints);
  TGraph gAZ(NPoints);

  TGraph2D g2D(NPoints);
  g2D.SetTitle("Electron Trajectory");
  g2D.GetXaxis()->SetTitle("X");
  g2D.GetYaxis()->SetTitle("Y");
  g2D.GetZaxis()->SetTitle("Z");

  std::vector<TVector3D> X, V, A;


  for (int i = 0; i != NPoints; ++i) {
    if (i % 1000000 == 0) {
      std::cout << "Step: " << i << std::endl;
    }
    double t = h * i;

    derivs(0, x, dxdt);

    X.push_back(TVector3D(x[0], x[2], x[4]));
    V.push_back(TVector3D(x[1], x[3], x[5]));
    A.push_back(TVector3D(dxdt[1], dxdt[3], dxdt[5]));

    gXZ.SetPoint(i, x[4], x[0]);
    gYZ.SetPoint(i, x[4], x[2]);
    gYX.SetPoint(i, x[0], x[2]);

    gVX.SetPoint(i, kC_SI * t, x[1] / kC_SI);
    gVY.SetPoint(i, kC_SI * t, x[3] / kC_SI);
    gVZ.SetPoint(i, kC_SI * t, x[5] / kC_SI);
    gAX.SetPoint(i, kC_SI * t, dxdt[1]);
    gAY.SetPoint(i, kC_SI * t, dxdt[3]);
    gAZ.SetPoint(i, kC_SI * t, dxdt[5]);

    g2D.SetPoint(i, x[0], x[2], x[4]);

    rk4(x, dxdt, N, t, h, x, derivs);
  }

  std::cout << x[4] << std::endl;
  std::cout << kC_SI * h * (NPoints-1) << std::endl;

  TCanvas c;
  c.cd();

  gXZ.Draw("AP");
  c.SaveAs("gXZ.png");

  gYZ.Draw("AP");
  c.SaveAs("gYZ.png");

  TFrame* Frame = (TFrame*) c.DrawFrame(-10e-6, -10e-6, 10e-6, 10e-6);
  Frame->Draw();
  gYX.Draw("P");
  c.SaveAs("gYX.png");


  c.cd();
  gVX.Draw("AP");
  c.SaveAs("gVX.png");
  gVY.Draw("AP");
  c.SaveAs("gVY.png");
  gVZ.Draw("AP");
  c.SaveAs("gVZ.png");

  gAX.Draw("AP");
  c.SaveAs("gAX.png");
  gAY.Draw("AP");
  c.SaveAs("gAY.png");
  gAZ.Draw("AP");
  c.SaveAs("gAZ.png");

  g2D.Draw("APL");
  c.SaveAs("g2D.png");



  return 0;
}



int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InBFieldFile]" << std::endl;
    return 1;
  }


  TBF = new TBField3DZRegularized(argv[1]);
  TBF->SaveAs("del2.dat");



  RK4Test();

  delete TBF;

  return 0;
}

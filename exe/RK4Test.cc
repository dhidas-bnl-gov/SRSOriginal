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
TVector3D GetValueAtTime (std::vector<TVector3D> const&, double const, double const, double const);
double RetardedTime (TVector3D const&, TVector3D const&, double const);
TVector3D ElectricField (TVector3D const& Observer, std::vector<TVector3D> const& X, std::vector<TVector3D> const& V, std::vector<TVector3D> const& A,  double const Time);


TBField3DZRegularized* TBF;

TVector3D GetValueAtTime (std::vector<TVector3D> const& V, double const T, double const T0, double const TStep)
{
  double const TD = T - T0;
  int const FirstBin = TD / TStep;

  if (FirstBin - 1 > (int) V.size()) {
    std::cout << "shit" << std::endl;
    throw;
  }


  TVector3D const A = (V[FirstBin + 1] - V[FirstBin]);

  double const Frac = (TD / TStep - ((int) (TD / TStep)));



  return V[FirstBin] + (A * Frac);
}


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
double const kEpsilon0 = 8.854187817E-12;
double const kMu0      = 1.2566370614E-6;

double const BetaX = 0;
double const BetaY = 0;



double const kC_SI = 299792458.;
//double const XStart = -1.8076921081543;
double const XStart = -1.5;
double const XStop  =  1.5;



double RetardedTime (TVector3D const& Object, TVector3D const& Observer, double const Time)
{
  return Time - (Object - Observer).Mag() / kC_SI;
}





double FutureTime (TVector3D const& Object, TVector3D const& Observer, double const RetardedTime)
{
  return RetardedTime + (Object - Observer).Mag() / kC_SI;
}




TVector3D ElectricField (TVector3D const& Observer, std::vector<TVector3D> const& X, std::vector<TVector3D> const& V, std::vector<TVector3D> const& A,  double const RTime, double const RT0, double const RTStep)
{
  TVector3D const XP = GetValueAtTime(X, RTime, RT0, RTStep);
  TVector3D const VP = GetValueAtTime(V, RTime, RT0, RTStep);
  TVector3D const AP = GetValueAtTime(A, RTime, RT0, RTStep);


  TVector3D const N = (Observer - XP).UnitVector();
  TVector3D const R = Observer - XP;
  TVector3D const B = VP / kC_SI;
  TVector3D const Bp = AP / kC_SI;

  double    const      Mult = kECharge / (4.0 * kPI * kEpsilon0) * pow(1.0 / (1.0 - N.Dot(VP) / kC_SI), 3);
  TVector3D const NearField = ((1.0 - VP.Mag2() / (kC_SI * kC_SI)) * (N - VP / kC_SI)) / (Observer - XP).Mag2();
  TVector3D const  FarField = (1.0/kC_SI) * (N.Cross(N - VP/kC_SI).Cross(AP / kC_SI));

  return TVector3D(Mult * (NearField + FarField));
}







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


  int const NPointsForward = 20000;
  int const NPointsBack = 0.5 / ((XStop - XStart) / (NPointsForward - 1));
  double const h = (XStop - XStart) / (BetaZ * kC_SI) / (NPointsForward - 1);

  int const NPointsTotal = NPointsForward + NPointsBack;


  std::vector<TVector3D> X, V, A;


  for (int i = 0; i != NPointsForward; ++i) {
    if (i % 1000000 == 0) {
      std::cout << "Step: " << i << std::endl;
    }
    double t = h * i;

    derivs(0, x, dxdt);

    X.push_back(TVector3D(x[0], x[2], x[4]));
    V.push_back(TVector3D(x[1], x[3], x[5]));
    A.push_back(TVector3D(dxdt[1], dxdt[3], dxdt[5]));


    //g2D.SetPoint(i, x[0], x[2], x[4]);

    rk4(x, dxdt, N, t, h, x, derivs);
  }





  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  x[3] = 0.0;
  x[4] = XStart;
  x[5] = -BetaZ * kC_SI;

  std::reverse(X.begin(), X.end());
  std::reverse(V.begin(), V.end());
  std::reverse(A.begin(), A.end());


  for (int i = 0; i != NPointsBack; ++i) {
    if (i % 1000 == 0) {
      std::cout << "StepBack: " << i << std::endl;
    }
    double t = h * i;

    derivs(0, x, dxdt);

    X.push_back(TVector3D(x[0], x[2], x[4]));
    V.push_back(TVector3D(-x[1], -x[3], -x[5]));
    A.push_back(TVector3D(-dxdt[1], -dxdt[3], -dxdt[5]));


    rk4(x, dxdt, N, t, h, x, derivs);
  }

  std::reverse(X.begin(), X.end());
  std::reverse(V.begin(), V.end());
  std::reverse(A.begin(), A.end());







  // Graphs for viewing trajectory
  TGraph gXZ(NPointsTotal);
  TGraph gYZ(NPointsTotal);
  TGraph gYX(NPointsTotal);

  TGraph gVX(NPointsTotal);
  TGraph gVY(NPointsTotal);
  TGraph gVZ(NPointsTotal);

  TGraph gAX(NPointsTotal);
  TGraph gAY(NPointsTotal);
  TGraph gAZ(NPointsTotal);

  TGraph2D g2D(NPointsTotal);
  g2D.SetTitle("Electron Trajectory");
  g2D.GetXaxis()->SetTitle("X");
  g2D.GetYaxis()->SetTitle("Y");
  g2D.GetZaxis()->SetTitle("Z");



  double const StartTime = h * ((double) (-NPointsBack));
  std::cout << "EndTime: " << h * ((double) (NPointsForward)) << std::endl;


  for (size_t i = 0; i != X.size(); ++i) {
    double const t = h * ((double) ((int) i - NPointsBack));
    gXZ.SetPoint(i, X[i].GetZ(), X[i].GetX());
    gYZ.SetPoint(i, X[i].GetZ(), X[i].GetY());
    gYX.SetPoint(i, X[i].GetX(), X[i].GetY());

    gVX.SetPoint(i, kC_SI * t, V[i].GetX() / kC_SI);
    gVY.SetPoint(i, kC_SI * t, V[i].GetY() / kC_SI);
    gVZ.SetPoint(i, kC_SI * t, V[i].GetZ() / kC_SI);
    gAX.SetPoint(i, kC_SI * t, A[i].GetX());
    gAY.SetPoint(i, kC_SI * t, A[i].GetY());
    gAZ.SetPoint(i, kC_SI * t, A[i].GetZ());

    g2D.SetPoint(i, X[i].GetX(), X[i].GetY(), X[i].GetZ());
  }

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




  TVector3D O(0, 0, 10);
  TVector3D const X0(X[0]);
  double const TT = h * ((double) (-NPointsBack));

  double const TFirstLight = (O - X0).Mag() / kC_SI;
  std::cout << "StartTime   " << StartTime << std::endl;
  std::cout << "TFirstLight " << TFirstLight << std::endl;
  std::cout << "RetardedTime FirstLight " << RetardedTime(X0, O, TFirstLight) << std::endl;

  TVector3D XT = GetValueAtTime(X, StartTime + h*20000.5, StartTime, h);
  std::cout << X[20000].GetZ() << std::endl;
  std::cout << XT.GetZ() << std::endl;

  /*
  for (double t = h*NPointsBack*0.9; t < (NPointsForward-1)*h ; t += h) {
    TVector3D const Ele = ElectricField(O, X, V, A, t, StartTime, h);
    //printf("%40.30E    %40.30E\n", t + (O - GetValueAtTime(X, t, StartTime, h)).Mag() / kC_SI , pow(Ele.GetX(), 2) + pow(Ele.GetY(), 2));
    printf("%40.30E    %40.30E\n", t + (O - GetValueAtTime(X, t, StartTime, h)).Mag() / kC_SI , Ele.GetX());
  }
  */

  O.SetZ(50);
  O.SetX(-0.05);
  TVector3D const XStep(0.001, 0, 0);
  TVector3D const YStep(0, 0.001, 0);
  for (int ix = 0; ix != 100; ++ix) {
    O.SetY(-0.05);
    O += XStep;
    for (int iy = 0; iy != 100; ++iy) {
      O += YStep;
      double Sum = 0.0;
      int Count = 0;
      for (double t = h*NPointsBack*0.9; t < (NPointsForward-1)*h ; t += h) {
        TVector3D const Ele = ElectricField(O, X, V, A, t, StartTime, h);
        Sum += pow(Ele.GetX(), 2) + pow(Ele.GetY(), 2);
        ++Count;
      }
      printf("%20.15E %20.15E    %40.30E\n", O.GetX(), O.GetY(),  Sum / (kMu0 / h));
    }
  }

  /*
  O.SetXYZ(0, -0.05, 50);
  TVector3D const XStep(0, 0.001, 0);
  for (int ix = 0; ix != 100; ++ix) {
    O += XStep;
    double Sum = 0.0;
    for (double t = h*NPointsBack*0.9; t < (NPointsForward-1)*h ; t += h) {
      TVector3D const Ele = ElectricField(O, X, V, A, t, StartTime, h);
      Sum += pow(Ele.GetX(), 2) + pow(Ele.GetY(), 2);
    }
    printf("%20.15E    %40.30E\n", O.GetX(), Sum);
  }
  */
  


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

////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Feb 10 11:14:09 EST 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "TBFieldUniformB.h"
#include "TBField3DZRegularized.h"
#include "TBFieldSquareWave.h"
#include "TBFieldIdeal1D.h"
#include "TVector3D.h"

#include "TGraph.h"
#include "TGraph2D.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFrame.h"





double const EEnergy = 3.0;
double const EMass   = 0.000511;
double const EGamma  = EEnergy / EMass;
double const EBeta   = sqrt(1.0 - 1.0/(EGamma*EGamma));

double const kPI = 3.14159265359;
double const k2PI = 2 * 3.14159265359;
double const BetaZ = EBeta;
double const Beta = BetaZ;
double const Gamma = EGamma;  //1.0 / sqrt(1.0 - Beta*Beta);
double const kEMass = 9.10938356E-31;
double const kECharge = -1.60217662E-19;
double const kEpsilon0 = 8.854187817E-12;
double const kMu0      = 1.2566370614E-6;

double const BetaX = 0;
double const BetaY = 0;



double const kC_SI = 299792458.;
//double const XStart = -1.8076921081543;
double const XStart = -2.5;
double const XStop  =  2.5;




int RK4Test ();
void rk4(double y[], double dydx[], int n, double x, double h, double yout[], void (*derivs)(double, double [], double []));
void derivs(double t, double y[], double dydt[]);
TVector3D GetValueAtTime (std::vector<TVector3D> const&, double const, double const, double const);
double RetardedTime (TVector3D const&, TVector3D const&, double const);
TVector3D ElectricField (TVector3D const& Observer, std::vector<TVector3D> const& X, std::vector<TVector3D> const& V, std::vector<TVector3D> const& A,  double const Time);

TVector3D PoyntingVector(TVector3D const&, TVector3D const&, TVector3D const&);


TBField* TBF;





TVector3D PoyntingVector(TVector3D const& El, TVector3D const& Particle, TVector3D const& Observer)
{
  TVector3D const N = (Observer - Particle).UnitVector();

  return 1. / (kMu0 * kC_SI) * (El.Mag2() * N - ( (El.Dot(N) * El) ));
}





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

  double    const      Mult = kECharge / (4.0 * kPI * kEpsilon0) * pow(1.0 / (1.0 - N.Dot(VP) / kC_SI), 3);
  //TVector3D const NearField = ((1.0 - VP.Mag2() / (kC_SI * kC_SI)) * (N - VP / kC_SI)) / (Observer - XP).Mag2();

  TVector3D const  FarField = (1.0 / kC_SI) * (N.Cross(  (N - VP/kC_SI).Cross(AP / kC_SI))  ) / (Observer - XP).Mag();

  //return TVector3D(Mult * (NearField + FarField));
  return TVector3D(Mult * FarField);
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


















int RK4Test ()
{


  int const N = 6;
  double x[N] = {0.0, 0.0, 0.0, 0.0, (double) XStart, (double) (BetaZ * kC_SI)};
  double dxdt[N];


  int const NPointsForward = 45001;
  int const NPointsBack = 0; //0.5 / ((XStop - XStart) / (NPointsForward - 1));
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
  std::cout << "StartTime: " << StartTime << std::endl;
  std::cout << "EndTime:   " << h * ((double) (NPointsForward)) << std::endl;


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




  TVector3D O(0.0005, 0.0000, 30);

  if (true) {
    TBField3DZ EleField;

    FILE* of1 = fopen("out1.dat", "w");
    //of1 << std::scientific;

    //for (double t = -NPointsBack * h; t < (NPointsForward-1)*h ; t += h) {
    for (double t = 0; t < (NPointsForward-1)*h ; t += h) {

      TVector3D const Ele = ElectricField(O, X, V, A, t, StartTime, h);
      EleField.Add(t + (O - GetValueAtTime(X, t, StartTime, h)).Mag() / kC_SI, Ele.GetX(), Ele.GetY(), Ele.GetZ());
      //printf("%40.30E    %40.30E\n", t + (O - GetValueAtTime(X, t, StartTime, h)).Mag() / kC_SI , pow(Ele.GetX(), 2) + pow(Ele.GetY(), 2));
    }

    TBField3DZRegularized EleFieldReg(EleField);



    //fprintf(of1, "%40.30E    %40.30E\n", t + (O - GetValueAtTime(X, t, StartTime, h)).Mag() / kC_SI , Ele.GetZ());
    fclose(of1);
  }





  if (false) {
    std::ofstream of2("out2.dat");
    of2 << std::scientific;

    int const npx = 40;
    double const xstart = -0.02;
    double const xstop = 0.02;
    double const xstep = (xstop-xstart) / npx;
    TVector3D const XStep(xstep, 0, 0);

    O.SetXYZ(xstart, 0, 30);

    double const tp0 = (O - GetValueAtTime(X, 0, StartTime, h)).Mag() / kC_SI;
    double const tp1 = (NPointsForward-1)*h + (O - GetValueAtTime(X, (NPointsForward-1)*h, StartTime, h)).Mag() / kC_SI;
    printf("%40.30E\n", tp1 - tp0);
    for (int ix = 0; ix != npx; ++ix) {
      O += XStep;

      double Sum = 0.0;
      int Count = 0;
      for (double t = 0; t < (NPointsForward-1)*h ; t += h) {
        TVector3D const Ele = ElectricField(O, X, V, A, t, StartTime, h);
        TVector3D const Particle = GetValueAtTime(X, t, StartTime, h);
        Sum += PoyntingVector(Ele, Particle, O).GetZ(); // Same as Mag() but faster in this case
        ++Count;
      }
      double const hp = (tp1 - tp0) / (double) Count;
      Sum *= hp;  // Gives time average in lab (assuming burst is < 1s)
      //Sum *= h;  // Gives time average (assuming burst is < 1s)
      Sum *= (0.5 / -kECharge); // from single electron per second to 500 mA
      Sum *= 1e-6; // from W/m^2 to W/mm^2
      //Sum /= (kMu0 * kC_SI);
      of2 << O.GetX() << " " << Sum << "\n";
    }
    of2.close();
  }


  if (false) {
    std::ofstream of3("out3.dat");
    of3 << std::scientific;

    int const npx = 40;
    double const xstart = -0.03;
    double const xstop = 0.03;
    double const xstep = (xstop-xstart) / npx;
    TVector3D const XStep(xstep, 0, 0);

    int const npy = 1;
    double const ystart = 0;
    double const ystop = 0.03;
    double const ystep = (ystop-ystart) / npy;
    TVector3D const YStep(0, ystep, 0);
    O.SetXYZ(xstart, ystart, 30);

    double TotalSum = 0.0;

    double const tp0 = (O - GetValueAtTime(X, 0, StartTime, h)).Mag() / kC_SI;
    double const tp1 = (NPointsForward-1)*h + (O - GetValueAtTime(X, (NPointsForward-1)*h, StartTime, h)).Mag() / kC_SI;
    for (int ix = 0; ix != npx; ++ix) {
      O += XStep;

      O.SetY(ystart);
      for (int iy = 0; iy != npy; ++iy) {
        O += YStep;

        double Sum = 0.0;
        int Count = 0;
        for (double t = 0; t < (NPointsForward-1)*h ; t += h) {
          TVector3D const Ele = ElectricField(O, X, V, A, t, StartTime, h);
          TVector3D const Particle = GetValueAtTime(X, t, StartTime, h);
          Sum += PoyntingVector(Ele, Particle, O).GetZ();
          ++Count;
        }
        double const hp = (tp1 - tp0) / (double) Count;
        Sum *= hp;  // Gives time average in lab (assuming burst is < 1s)
        //Sum *= h;  // Gives time average (assuming burst is < 1s)
        Sum *= (0.5 / -kECharge); // from single electron per second to 500 mA
        Sum *= 1e-6; // from W/m^2 to W/mm^2
        of3 << O.GetX() << " " << O.GetY() << " " << Sum << "\n";

        TotalSum += Sum;
      }
    }
    of3.close();

    std::cout << "Power: " << TotalSum * 1e6 * xstep * ystep << std::endl;
  }


  return 0;
}



int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InBFieldFile]" << std::endl;
    return 1;
  }


  //TBF =(TBField*) new TBField3DZRegularized(argv[1]);
  //TBF =(TBField*) new TBFieldSquareWave(0.200, 11, 0, 1.0);
  TBF =(TBField*) new TBFieldIdeal1D(0.057, 22, 0, 0.836);
  //TBF =(TBField*) new TBFieldUniformB(0, 0, 0);



  RK4Test();

  delete TBF;

  return 0;
}

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
#include <complex>

#include "TSurfacePoints_RectangleSimple.h"
#include "TSurfacePoints_BBoxSimple.h"
#include "TBFieldUniformB.h"
#include "TBField3DZRegularized.h"
#include "TBFieldSquareWave.h"
#include "TBFieldIdeal1D.h"
#include "TVector3D.h"
#include "TVector3DC.h"

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
double const kh = 6.62607004e-34;

double const BetaX = 0;
double const BetaY = 0;




double const kC_SI = 299792458.;
//double const XStart = -1.8076921081543;
double const XStart = -1.5;
double const XStop  =  1.5;




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






double PowerDensityIntegrand (TVector3D const& X, TVector3D const& V, TVector3D const& A, TVector3D const& O, TVector3D const& U)
{
  TVector3D const N = (O - X).UnitVector();

  double const Numerator = N.Cross( ( (N - (V / kC_SI)).Cross((A / kC_SI)) ) ).Dot(U);
  double const Denominator = pow(1 - (V / kC_SI).Dot(N), 5);


  return Numerator*Numerator / Denominator / (O - X).Mag2();
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


  int const NPointsForward = 20001;
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







  // Beam current in Amps
  double const Current = 0.500;


  if (false) {
    std::ofstream oB("out_ZBy.dat");
    oB << std::scientific;
    for (double z = XStart; z <= XStop; z += 0.0001) {
      //oB << 0 << "\t" << TBF->GetBy(0, 0, z) << "\t" << 0 << std::endl;
      oB << z << "\t" << TBF->GetBy(0, 0, z) << std::endl;
    }
    oB.close();
    exit(0);
  }







  if (false) {
    std::ofstream ofFlux("out_flux.dat");
    ofFlux << std::scientific;

    TVector3D Obs(0, 0.000, 300);
    TSurfacePoints_RectangleSimple Surface("XY", 51, 51, 160e-3, 30e-3, 0, 0, 30, 1);


    double const C0 = kECharge / (4 * kPI * kC_SI * kEpsilon0);
    std::complex<double> const I(0, 1);

    double const Energy_eV = 5000;//188;//135;//188; //565;//188;
    double const iw = Energy_eV * k2PI / 4.1357e-15;

    std::complex<double> const C1(0, C0 * iw);

    for (size_t i = 0; i != Surface.GetNPoints(); ++i) {
      if (i % 100 == 0) {
        std::cout << ( (int) ((double) i / (double) Surface.GetNPoints() * 100.) ) << " \% done" << std::endl;
      }

      TVector3D Obs = Surface.GetPoint(i).GetPoint();
      TVector3DC Normal = Surface.GetPoint(i).GetNormal();


      TVector3DC SumE(0, 0, 0);
      TVector3DC SumB(0, 0, 0);

      for (int i = 0; i != NPointsForward; ++i) {
        TVector3D const R = Obs - X[i];
        TVector3D const N = R.UnitVector();
        double const D = R.Mag();
        std::complex<double> Exponent(0, iw * (h * i + D / kC_SI));

        TVector3DC const ThisEw = (TVector3DC(V[i]) / kC_SI - TVector3DC(N) * ( std::complex<double>(1, 0) + (I * kC_SI / (iw * D)))) / D * std::exp(Exponent) * h;
        TVector3DC const ThisBw = ((TVector3DC(V[i]) / kC_SI).Cross( TVector3DC(N) * ( std::complex<double>(1, 0) + (I * kC_SI / (iw * D)))) / D * std::exp(Exponent) * h );

        SumE += ThisEw;
        SumB -= ThisBw / kC_SI;

      }


      SumE *= C1;
      SumB *= C1;

      ofFlux << Obs.GetX() << " " << Obs.GetY() << " " <<  8 * kPI * kPI * kEpsilon0 * kC_SI * kC_SI * Current / (kh * fabs(kECharge)) * SumE.Cross(SumB.CC()).Dot( Normal ).real() * 1e-6 << std::endl; // fux per m,^2
    }

    ofFlux.close();
    exit(0);

  }












  if (true) {
    // Spec according to Hoffman
    std::ofstream ofSpec("out_spec.dat");
    ofSpec << std::scientific;

    TVector3D Obs(0, 0.000, 190);



    double const C0 = kECharge / (4 * kPI * kC_SI * kEpsilon0 * sqrt(k2PI));
    std::complex<double> const I(0, 1);

    int const NEPoints = 1000;
    double const EStart = 100;
    double const EStop  = 2000;

    double const wStart = EStart * k2PI / 4.1357e-15;
    double const wStop  = EStop * k2PI / 4.1357e-15;

    double const wStepSize = (wStop - wStart) / (NEPoints - 1);

    for (double iw = wStart; iw <= wStop; iw += wStepSize) {



      std::complex<double> SumX(0, 0);
      std::complex<double> SumY(0, 0);
      std::complex<double> SumZ(0, 0);

      TVector3DC SumE(0, 0, 0);

      for (int i = 0; i != NPointsForward; ++i) {
        TVector3D const R = Obs - X[i];
        TVector3D const N = R.UnitVector();
        double const D = R.Mag();
        std::complex<double> Exponent(0, -iw * (h * i + D / kC_SI));

        TVector3DC const ThisEw = ( N.Cross( (N - V[i] / kC_SI).Cross(A[i] / kC_SI) ) ) / ( D * pow(1 - N.Dot(V[i] / kC_SI), 2) ) * std::exp(Exponent) * h;

        SumE += ThisEw;

      }


      SumE *= C0;

      //ofSpec << iw * 4.1357e-15 / k2PI << "  " <<  2 * k2PI * Current / (kh * fabs(kECharge) * kMu0 * kC_SI) *  SumE.Dot( SumE.CC() ).real() * pow(Obs.GetZ(), 2) << std::endl;
      ofSpec << iw * 4.1357e-15 / k2PI << "  " <<  2 * k2PI * Current / (kh * fabs(kECharge) * kMu0 * kC_SI) *  SumE.Dot( SumE.CC() ).real()  * 1e-6  * 0.001 << std::endl;

    }

    ofSpec.close();

    exit(0);
  }






















  if (false) {
    // Spec according to Und, Wigg, and Apps
    std::ofstream ofSpec("out_spec.dat");
    ofSpec << std::scientific;

    TVector3D Obs(0, 0.000, 30);



    double const C0 = kECharge / (4 * kPI * kC_SI * kEpsilon0);
    std::complex<double> const I(0, 1);

    int const NEPoints = 1000;
    double const EStart = 100;
    double const EStop  = 2000;

    double const wStart = EStart * k2PI / 4.1357e-15;
    double const wStop  = EStop * k2PI / 4.1357e-15;

    double const wStepSize = (wStop - wStart) / (NEPoints - 1);

    for (double iw = wStart; iw <= wStop; iw += wStepSize) {
      std::complex<double> const C1(0, C0 * iw);



      std::complex<double> SumX(0, 0);
      std::complex<double> SumY(0, 0);
      std::complex<double> SumZ(0, 0);

      TVector3DC SumE(0, 0, 0);
      TVector3DC SumB(0, 0, 0);

      for (int i = 0; i != NPointsForward; ++i) {
        TVector3D const R = Obs - X[i];
        TVector3D const N = R.UnitVector();
        double const D = R.Mag();
        std::complex<double> Exponent(0, iw * (h * i + D / kC_SI));

        TVector3DC const ThisEw = (TVector3DC(V[i]) / kC_SI - (TVector3DC(N) * ( std::complex<double>(1, 0) + (I * kC_SI / (iw * D))))) / D * std::exp(Exponent) * h;
        TVector3DC const ThisBw = ((TVector3DC(V[i]) / kC_SI).Cross( TVector3DC(N) * ( std::complex<double>(1, 0) + (I * kC_SI / (iw * D)))) / D * std::exp(Exponent) * h );

        SumE += ThisEw;
        SumB += ThisBw;

      }


      SumE *= C1;
      SumB *= -C1 / kC_SI;

      ofSpec << iw * 4.1357e-15 / k2PI << "  " <<  8 * kPI * kPI * kEpsilon0 * kC_SI * kC_SI * Current / (kh * fabs(kECharge)) * ( SumE.Cross(SumB.CC()).Dot( TVector3DC(0, 0, 1) ) ).real() * 1e-6 << std::endl;




    }

    ofSpec.close();

    exit(0);
  }



  if (false) {

    TVector3D Obs(0, 0.000, 30);

    double const ConvmRad = pow(Obs.GetZ() * 0.001, 2);


    double const C0 = kECharge / (4 * kPI * kC_SI * kEpsilon0);
    std::complex<double> const I(0, 1);

    int const NEPoints = 10000;
    double const EStart = 1;
    double const EStop  = 50000;

    double const wStart = EStart * k2PI / 4.1357e-15;
    double const wStop  = EStop * k2PI / 4.1357e-15;

    double const wStepSize = (wStop - wStart) / (NEPoints - 1);

    double TotalPD = 0;

    for (double iw = wStart; iw <= wStop; iw += wStepSize) {
      std::complex<double> const C1(0, C0 * iw);



      TVector3DC SumE(0, 0, 0);
      TVector3DC SumB(0, 0, 0);

      for (int i = 0; i != NPointsForward; ++i) {
        TVector3D const R = Obs - X[i];
        TVector3D const N = R.UnitVector();
        double const D = R.Mag();
        std::complex<double> Exponent(0, iw * (h * i + D / kC_SI));


        TVector3DC const ThisEw = (TVector3DC(V[i]) / kC_SI - TVector3DC(N) * ( std::complex<double>(1, 0) + (I * kC_SI / (iw * D)))) / D * std::exp(Exponent) * h;
        TVector3DC const ThisBw = ((TVector3DC(V[i]) / kC_SI).Cross( TVector3DC(N) * ( std::complex<double>(1, 0) + (I * kC_SI / (iw * D)))) / D * std::exp(Exponent) * h );

        SumE += ThisEw;
        SumB -= ThisBw / kC_SI;
      }


      SumE *= C1;
      SumB *= C1;

      TotalPD +=   4 * kPI * kEpsilon0 * kC_SI * kC_SI * Current / (fabs(kECharge)) * SumE.Cross(SumB.CC()).Dot( TVector3DC(0, 0, 1) ).real() * 1e-6 * wStepSize;



    }

    std::cout << "PD: " << TotalPD << std::endl;


    exit(0);
  }




  // Calculate total power out
  double TotalPower = 0;
  for (int i = 0; i != NPointsForward; ++i) {
    TotalPower += ( (A[i] / kC_SI).Mag2() - ( (V[i] / kC_SI).Cross(A[i] / kC_SI)).Mag2() ) * h;
  }
  TotalPower *= fabs(kECharge * Current) * pow(Gamma, 6) / (6 * kPI * kEpsilon0 * kC_SI);
  std::cout << "Total Power: " << TotalPower << std::endl;

  TotalPower = 0;
  for (int i = 0; i != NPointsForward; ++i) {
    TotalPower += ( (A[i] / kC_SI).Mag2() - ( (V[i] / kC_SI).Cross(A[i] / kC_SI)).Mag2() ) * h;
  }
  TotalPower *= fabs(kECharge * Current) * pow(Gamma, 6) / (6 * kPI * kEpsilon0 * kC_SI);
  std::cout << "Total Power: " << TotalPower << std::endl;






  if (false) {
    std::vector<std::ofstream> of(6);
    of[0].open("BBoxSimple_PXY.dat");
    of[1].open("BBoxSimple_MXY.dat");
    of[2].open("BBoxSimple_PXZ.dat");
    of[3].open("BBoxSimple_MXZ.dat");
    of[4].open("BBoxSimple_PYZ.dat");
    of[5].open("BBoxSimple_MYZ.dat");
    for (int i = 0; i != 6; ++i) {
      of[i] << std::scientific;
    }
    TSurfacePoints_BBoxSimple Surface(0, 0, 0, 0.060, 0.060, 60, 50, 50, 1);

    for (size_t io = 0; io < Surface.GetNPoints(); ++io) {
      TVector3D Obs = Surface.GetPoint(io).GetPoint();
      TVector3D Normal = Surface.GetPoint(io).GetNormal();

      double Sum = 0;

      for (int i = 0; i != NPointsForward ; ++i) {
        TVector3D const N1 = (Obs - X[i]).UnitVector();
        TVector3D const N2 = N1.Cross(TVector3D(1, 0, 0)).UnitVector();
        TVector3D const N3 = N1.Cross(N2).UnitVector();

        Sum += PowerDensityIntegrand(X[i], V[i], A[i], Obs, N2) * N1.Dot(Normal);
        Sum += PowerDensityIntegrand(X[i], V[i], A[i], Obs, N3) * N1.Dot(Normal);
      }

      // Put into SI units
      Sum *= fabs(kECharge * Current) / (16 * kPI * kPI * kEpsilon0 * kC_SI) * h;

      Sum /= 1e6; // m^2 to mm^2


      of[Surface.GetSurfaceNumber(io)] << Surface.GetX1(io) << " " << Surface.GetX2(io) << " " << Sum << "\n";

    }

    for (std::vector<std::ofstream>::iterator it = of.begin(); it != of.end(); ++it) {
      it->close();
    }
    of.clear();
    exit(0);
  }




  double TotalSum = 0;
  //TSurfacePoints_RectangleSimple Surface("XY", 101, 101, 0.16, 0.026, 0, 0, 30, 1);
  TSurfacePoints_RectangleSimple Surface("XY", 101, 101, 0.06, 0.06, 0, 0, 30, 1);
  //TSurfacePoints_RectangleSimple Surface("XZ", 51, 51, 0.04, 3, 0, 0.004, 0, 1);
  std::ofstream ofS("out.dat");
  ofS << std::scientific;
  for (size_t io = 0; io != Surface.GetNPoints(); ++io) {

    TVector3D Obs = Surface.GetPoint(io).GetPoint();
    TVector3D Normal = Surface.GetPoint(io).GetNormal();

    double Sum = 0;

    for (int i = 0; i != NPointsForward ; ++i) {
      TVector3D const N1 = (Obs - X[i]).UnitVector();
      TVector3D const N2 = N1.Cross(TVector3D(1, 0, 0)).UnitVector();
      TVector3D const N3 = N1.Cross(N2).UnitVector();

      Sum += PowerDensityIntegrand(X[i], V[i], A[i], Obs, N2) * N1.Dot(Normal);
      Sum += PowerDensityIntegrand(X[i], V[i], A[i], Obs, N3) * N1.Dot(Normal);
    }

    // Put into SI units
    Sum *= fabs(kECharge * Current) / (16 * kPI * kPI * kEpsilon0 * kC_SI) * h;

    Sum /= 1e6; // m^2 to mm^2


    ofS << Surface.GetX1(io) << " " << Surface.GetX2(io) << " " << Sum << "\n";

    TotalSum += Sum;

  }
  std::cout << "Power: " << TotalSum * Surface.GetElementArea() * 1e6 << std::endl;
  ofS.close();
  exit(0);

  if (false) {
    std::ofstream of2("out2.dat");
    of2 << std::scientific;

    TVector3D O(0.0000, 0.0000, 30);

    int const npx = 40;
    double const xstart = -0.03;
    double const xstop  =  0.03;
    double const xstep  = (xstop-xstart) / (npx - 1);
    TVector3D const XStep(xstep, 0, 0);

    O.SetXYZ(xstart, 0, 30);

    for (int ix = 0; ix != npx; ++ix) {
      //TVector3D const U = (O - TVector3D(0, 0, 0)).UnitVector();
      TVector3D const U = TVector3D(1, 0, 0);

      double Sum = 0.0;
      for (int i = 0; i != NPointsForward ; ++i) {
        Sum += PowerDensityIntegrand(X[i], V[i], A[i], O, U);
      }
      // Put into SI units
      Sum *= fabs(kECharge * Current) / (16 * kPI * kPI * kEpsilon0 * kC_SI) * h;

      Sum *= 1e3; // W to mW
      Sum /= 1e6; // m^2 to mm^2

      of2 << O.GetX() << "  " << Sum << std::endl;
      O += XStep;
    }

    of2.close();
  }


  if (false) {
    std::ofstream of3("out_plusZ.dat");
    of3 << std::scientific;

    TVector3D O(0.0000, 0.0000, 30);

    int const npx = 101;
    double const xstart = -0.02;
    double const xstop  =  0.02;
    double const xstep  = (xstop-xstart) / (npx - 1);
    TVector3D const XStep(xstep, 0, 0);

    int const npy = 101;
    double const ystart = -0.004;
    double const ystop  =  0.004;
    double const ystep  = (ystop-ystart) / (npy - 1);
    TVector3D const YStep(0, ystep, 0);
    O.SetXYZ(xstart, ystart, 6);

    double TotalSum = 0.0;

    TVector3D const UX = TVector3D(1, 0, 0);
    TVector3D const UY = TVector3D(0, 1, 0);
    for (int ix = 0; ix != npx; ++ix) {

      O.SetY(ystart);
      for (int iy = 0; iy != npy; ++iy) {

        double Sum = 0;

        for (int i = 0; i != NPointsForward ; ++i) {
          TVector3D const N1 = (O - X[i]).UnitVector();
          TVector3D const N2 = N1.Cross(TVector3D(1, 0, 0)).UnitVector();
          TVector3D const N3 = N1.Cross(N2).UnitVector();

          Sum += PowerDensityIntegrand(X[i], V[i], A[i], O, N2) * N1.Dot(TVector3D(0, 0, 1));
          Sum += PowerDensityIntegrand(X[i], V[i], A[i], O, N3) * N1.Dot(TVector3D(0, 0, 1));
          //Sum += PowerDensityIntegrand(X[i], V[i], A[i], O, UX);
          //Sum += PowerDensityIntegrand(X[i], V[i], A[i], O, UY);
        }

        // Put into SI units
        Sum *= fabs(kECharge * Current) / (16 * kPI * kPI * kEpsilon0 * kC_SI) * h;

        //Sum *= 1e3; // W to mW
        Sum /= 1e6; // m^2 to mm^2



        of3 << O.GetX() << " " << O.GetY() << " " << Sum << "\n";

        TotalSum += Sum;
        O += YStep;
      }
      O += XStep;
    }
    of3.close();

    std::cout << "Power: " << TotalSum * xstep * ystep * 1e6 << std::endl;
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
  TBF =(TBField*) new TBFieldIdeal1D(0.03, 33, 0, 1.000);
  //TBF =(TBField*) new TBFieldUniformB(0, 0, 0);



  RK4Test();

  delete TBF;

  return 0;
}

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

#include "TSRS.h"
#include "TParticleTrajectoryPoints.h"
#include "TSurfacePoints_RectangleSimple.h"
#include "TSurfacePoints_BBoxSimple.h"
#include "TBFieldUniformB.h"
#include "TBField3DZRegularized.h"
#include "TBFieldSquareWave.h"
#include "TBFieldIdeal1D.h"
#include "TVector3D.h"
#include "TVector3DC.h"
#include "TSpectrumContainer.h"

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



TBField* TBF;













void CalculateSpectrumAtPoint (TParticleTrajectoryPoints const& T, TVector3D const& ObservationPoint, double const Current, TSpectrumContainer& S)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // T - Trajectory of particle
  // ObservationPoint - Observation Point
  // Current - beam current
  // S - TSpectrumContainer object which will be used to decide what energies to look at and save results to

  double const DeltaT = T.GetDeltaT();

  // Spec according to Hoffman

  size_t const NTPoints = T.GetNPoints();

  double const C0 = TSRS::Qe() / (TSRS::FourPi() * TSRS::C() * TSRS::Epsilon0() * TSRS::Sqrt2Pi());

  std::complex<double> const I(0, 1);

  size_t const NEPoints = S.GetNPoints();


  for (size_t i = 0; i != NEPoints; ++i) {
    double const Omega = S.GetAngularFrequency(i);

    TVector3DC SumE(0, 0, 0);

    for (int iT = 0; iT != NTPoints; ++iT) {
      TVector3D const& X = T.GetX(iT);
      TVector3D const& B = T.GetB(iT);
      TVector3D const& AoverC = T.GetAoverC(iT);

      TVector3D const R = ObservationPoint - X;
      TVector3D const N = R.UnitVector();
      double const D = R.Mag();



      std::complex<double> Exponent(0, -Omega * (DeltaT * iT + D / TSRS::C()));

      // This is in the far field calculation only
      // SumE = ( N.Cross( (N - B).Cross(AoverC) ) ) / ( D * pow(1 - N.Dot(B), 2) ) * std::exp(Exponent);

      // This is the full near + far field
      SumE += ( ( (1 - B.Mag2()) * (N - B) ) / ( D * D * (pow(1 - N.Dot(B), 2)) )
            + ( N.Cross( (N - B).Cross(AoverC) ) ) / ( D * pow(1 - N.Dot(B), 2) ) ) * std::exp(Exponent); // NF + FF

    }


    SumE *= C0 * DeltaT;

    S.SetFlux(i, TSRS::FourPi() * Current / (TSRS::H() * fabs(TSRS::Qe()) * TSRS::Mu0() * TSRS::C()) *  SumE.Dot( SumE.CC() ).real()  * 1e-6  * 0.001);

  }

  return;
}



void CalculateSpectrumAtPoint2 (TParticleTrajectoryPoints const& T, TVector3D const& ObservationPoint, double const Current, TSpectrumContainer& S)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // T - Trajectory of particle
  // ObservationPoint - Observation Point
  // Current - beam current
  // S - TSpectrumContainer object which will be used to decide what energies to look at and save results to

  double const DeltaT = T.GetDeltaT();

  // Spec according to Hoffman

  size_t const NTPoints = T.GetNPoints();

  double const C0 = TSRS::Qe() / (TSRS::FourPi() * TSRS::C() * TSRS::Epsilon0());

  std::complex<double> const I(0, 1);

  size_t const NEPoints = S.GetNPoints();

  for (size_t i = 0; i != NEPoints; ++i) {
    double const Omega = S.GetAngularFrequency(i);
    std::complex<double> const C1(0, C0 * Omega);

    TVector3DC SumE(0, 0, 0);
    TVector3DC SumB(0, 0, 0);

    for (int iT = 0; iT != NTPoints; ++iT) {
      TVector3D const& X = T.GetX(iT);
      TVector3D const& B = T.GetB(iT);

      TVector3D const R = ObservationPoint - X;
      TVector3D const N = R.UnitVector();
      double const D = R.Mag();


      std::complex<double> Exponent(0, Omega * (DeltaT * iT + D / TSRS::C()));

      TVector3DC const ThisEw = (TVector3DC(B) - (N * ( std::complex<double>(1, 0) + (I * TSRS::C() / (Omega * D))))) / D * std::exp(Exponent) * DeltaT;
      //TVector3DC const ThisBw = ((TVector3DC(B).Cross( TVector3DC(N)) * ( std::complex<double>(1, 0) + (I * TSRS::C() / (Omega * D)))) / D * std::exp(Exponent) * DeltaT );
      TVector3DC const ThisBw = ThisEw.Cross(N);


      SumE += ThisEw;
      SumB += ThisBw;

    }

    SumE *= C1;
    SumB *= -C1 / kC_SI;

    S.SetFlux(i, 8 * TSRS::Pi() * TSRS::Pi() * TSRS::Epsilon0() * TSRS::C() * TSRS::C() *  Current / (TSRS::H() * fabs(TSRS::Qe()) ) * ( SumE.Cross(SumB.CC()).Dot( TVector3DC(0, 0, 1) ) ).real() * 1e-6 * 0.001 );

  }

  return;
}






void CalculatePowerDensitySurface (TParticleTrajectoryPoints const& T, TSurfacePoints const& Surface, double const Current)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // T - Trajectory of particle
  // Surface - Observation Point
  // Current - beam current

  size_t const NTPoints = T.GetNPoints();
  double const DeltaT = T.GetDeltaT();

  TVector3D Numerator;
  double Denominator;

  std::ofstream of("out_pow.dat");
  of << std::scientific;

  for (size_t io = 0; io < Surface.GetNPoints(); ++io) {
    TVector3D Obs = Surface.GetPoint(io).GetPoint();
    TVector3D Normal = Surface.GetPoint(io).GetNormal();

    double Sum = 0;

    for (int iT = 0; iT != NTPoints ; ++iT) {
      TVector3D const& X = T.GetX(iT);
      TVector3D const& B = T.GetB(iT);
      TVector3D const& AoverC = T.GetAoverC(iT);


      TVector3D const N1 = (Obs - X).UnitVector();
      TVector3D const N2 = N1.Cross(TVector3D(1, 0, 0)).UnitVector();
      TVector3D const N3 = N1.Cross(N2).UnitVector();

      double const N1DotNormal = N1.Dot(Normal);

      Numerator = N1.Cross( ( (N1 - B).Cross((AoverC)) ) );
      Denominator = pow(1 - (B).Dot(N1), 5);

      Sum += pow(Numerator.Dot(N2), 2) / Denominator / (Obs - X).Mag2() * N1DotNormal;
      Sum += pow(Numerator.Dot(N3), 2) / Denominator / (Obs - X).Mag2() * N1DotNormal;
      //Sum += PowerDensityIntegrand(X[i], V[i], A[i], Obs, N2) * N1DotNormal;
      //Sum += PowerDensityIntegrand(X[i], V[i], A[i], Obs, N3) * N1DotNormal;
    }

    // Put into SI units
    Sum *= fabs(TSRS::Qe() * Current) / (16 * TSRS::Pi() * TSRS::Pi() * TSRS::Epsilon0() * TSRS::C()) * DeltaT;

    Sum /= 1e6; // m^2 to mm^2


    of << Surface.GetX1(io) << " " << Surface.GetX2(io) << " " << Sum << "\n";

  }
  of.close();
  exit(0);


  return;
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
  //dxdt[1] = -(kECharge * x[5]) / (Gamma * kEMass) * TBF->GetBy(x[0], x[2], x[4]) + (kECharge * x[3]) / (Gamma * kEMass) * TBF->GetBz(x[0], x[2], x[4]);
  dxdt[1] = TSRS::QeOverMe() / Gamma * (-x[5] * TBF->GetBy(x[0], x[2], x[4]) + x[3] * TBF->GetBz(x[0], x[2], x[4]));
  dxdt[2] = x[3];                                                                                                                            
  dxdt[3] = TSRS::QeOverMe() / Gamma * ( x[5] * TBF->GetBx(x[0], x[2], x[4]) - x[1] * TBF->GetBz(x[0], x[2], x[4]));
  dxdt[4] = x[5];                                                                                                                            
  dxdt[5] = TSRS::QeOverMe() / Gamma * ( x[1] * TBF->GetBy(x[0], x[2], x[4]) - x[3] * TBF->GetBx(x[0], x[2], x[4]));

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


  int const NPointsForward = 5001;
  int const NPointsBack = 0; //0.5 / ((XStop - XStart) / (NPointsForward - 1));
  double const h = (XStop - XStart) / (BetaZ * kC_SI) / (NPointsForward - 1);

  int const NPointsTotal = NPointsForward + NPointsBack;


  std::vector<TVector3D> X, V, A;
  TParticleTrajectoryPoints ParticleTrajectory(h);


  for (int i = 0; i != NPointsForward; ++i) {
    double t = h * i;

    derivs(0, x, dxdt);

    ParticleTrajectory.AddPoint(x[0], x[2], x[4], x[1] / TSRS::C(), x[3] / TSRS::C(), x[5] / TSRS::C(), dxdt[1] / TSRS::C(), dxdt[3] / TSRS::C(), dxdt[5] / TSRS::C());

    rk4(x, dxdt, N, t, h, x, derivs);
  }

  if (NPointsBack > 0) {

    // Reverse direction and prop back
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = XStart;
    x[5] = -BetaZ * kC_SI;

    ParticleTrajectory.ReverseArrays();

    for (int i = 0; i != NPointsBack; ++i) {
      if (i % 1000 == 0) {
        std::cout << "StepBack: " << i << std::endl;
      }
      double t = h * i;

      derivs(0, x, dxdt);

      ParticleTrajectory.AddPoint(x[0], x[2], x[4], -x[1] / TSRS::C(), -x[3] / TSRS::C(), -x[5] / TSRS::C(), -dxdt[1] / TSRS::C(), -dxdt[3] / TSRS::C(), -dxdt[5] / TSRS::C());

      rk4(x, dxdt, N, t, h, x, derivs);
    }

    ParticleTrajectory.ReverseArrays();

  }







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


  for (size_t i = 0; i != ParticleTrajectory.GetNPoints(); ++i) {
    double const t = h * ((double) ((int) i - NPointsBack));

    TVector3D X = ParticleTrajectory.GetX(i);
    TVector3D B = ParticleTrajectory.GetB(i);
    TVector3D A = ParticleTrajectory.GetA(i);


    gXZ.SetPoint(i, X.GetZ(), X.GetX());
    gYZ.SetPoint(i, X.GetZ(), X.GetY());
    gYX.SetPoint(i, X.GetX(), X.GetY());

    gVX.SetPoint(i, kC_SI * t, B.GetX());
    gVY.SetPoint(i, kC_SI * t, B.GetY());
    gVZ.SetPoint(i, kC_SI * t, B.GetZ());
    gAX.SetPoint(i, kC_SI * t, A.GetX());
    gAY.SetPoint(i, kC_SI * t, A.GetY());
    gAZ.SetPoint(i, kC_SI * t, A.GetZ());

    g2D.SetPoint(i, X.GetX(), X.GetY(), X.GetZ());
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







  TVector3D Obs(0.000, 0.000, 30);
  //TSpectrumContainer Spectrum(2000, 100, 2000);
  //CalculateSpectrumAtPoint2(ParticleTrajectory, Obs, 0.500, Spectrum);
  //CalculateSpectrumAtPoint(ParticleTrajectory, Obs, 0.500, Spectrum);
  //Spectrum.SaveToFile("out_spec.dat", "");

  // Calculate Power density on a surface
  TSurfacePoints_RectangleSimple Surface("XY", 101, 101, 0.06, 0.06, 0, 0, 30, 1);
  CalculatePowerDensitySurface(ParticleTrajectory, Surface, 0.500);


  // CalculateTotalPower()
  // CalculatePowerDensity()
  // CalculatePowerDensity() with gaussian convolution

  // CalculateFlux input SurfacePoint or whole Surface?
  //  How to keep data correlated...

  exit(0);





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







  if (true) {
    std::ofstream ofFlux("out_flux.dat");
    ofFlux << std::scientific;

    TSurfacePoints_RectangleSimple Surface("XY", 401, 401, 160e-3, 30e-3, 0, 0, 30, 1);

    double const C0 = kECharge / (4 * kPI * kC_SI * kEpsilon0 * sqrt(k2PI));
    std::complex<double> const I(0, 1);

    double const Energy_eV = 5000;//188;//135;//188; //565;//188;
    double const iw = Energy_eV * k2PI / 4.1357e-15;

    std::complex<double> const C1(0, C0 * iw);

    for (size_t ip = 0; ip != Surface.GetNPoints(); ++ip) {
      if (ip % 100 == 0) {
        std::cout << ( (int) ((double) ip / (double) Surface.GetNPoints() * 100.) ) << " \% done" << std::endl;
      }

      TVector3D Obs = Surface.GetPoint(ip).GetPoint();

      TVector3DC SumE(0, 0, 0);

      for (int i = 0; i != NPointsForward; ++i) {
        TVector3D const R = Obs - X[i];
        TVector3D const N = R.UnitVector();
        double const D = R.Mag();
        std::complex<double> Exponent(0, -iw * (h * i + D / kC_SI));

        // TVector3DC const ThisEw = ( N.Cross( (N - V[i] / kC_SI).Cross(A[i] / kC_SI) ) ) / ( D * pow(1 - N.Dot(V[i] / kC_SI), 2) ) * std::exp(Exponent) * h; // FF only
        TVector3DC const ThisEw = ( ( (1 - (V[i] / kC_SI).Mag2()) * (N - (V[i] / kC_SI)) ) / ( D * D * (pow(1 - N.Dot(V[i] / kC_SI), 2)) )
            + ( N.Cross( (N - V[i] / kC_SI).Cross(A[i] / kC_SI) ) ) / ( D * pow(1 - N.Dot(V[i] / kC_SI), 2) ) ) * std::exp(Exponent) * h; // NF + FF

        SumE += ThisEw;

      }


      SumE *= C0;

      //ofFlux << iw * 4.1357e-15 / k2PI << "  " <<  2 * k2PI * Current / (kh * fabs(kECharge) * kMu0 * kC_SI) *  SumE.Dot( SumE.CC() ).real()  * 1e-6  * 0.001 * Obs.Mag2() << std::endl;
      ofFlux << Surface.GetX1(ip) << "  " << Surface.GetX2(ip) << "  " <<  2 * k2PI * Current / (kh * fabs(kECharge) * kMu0 * kC_SI) *  SumE.Dot( SumE.CC() ).real()  * 1e-6  * 0.001 << std::endl;
    }

    ofFlux.close();
    exit(0);

  }












  if (false) {
    // Spec according to Hoffman
    std::ofstream ofSpec("out_spec.dat");
    ofSpec << std::scientific;

    TVector3D Obs(0.000, 0.000, 30);



    double const C0 = kECharge / (4 * kPI * kC_SI * kEpsilon0 * sqrt(k2PI));
    std::complex<double> const I(0, 1);

    int const NEPoints = 2000;
    double const EStart = 1000;
    double const EStop  = 10000;

    double const wStart = EStart * k2PI / 4.1357e-15;
    double const wStop  = EStop * k2PI / 4.1357e-15;

    double const wStepSize = (wStop - wStart) / (NEPoints - 1);

    for (double iw = wStart; iw <= wStop; iw += wStepSize) {




      TVector3DC SumE(0, 0, 0);

      for (int i = 0; i != NPointsForward; ++i) {
        TVector3D const R = Obs - X[i];
        TVector3D const N = R.UnitVector();
        double const D = R.Mag();
        std::complex<double> Exponent(0, -iw * (h * i + D / kC_SI));

        // TVector3DC const ThisEw = ( N.Cross( (N - V[i] / kC_SI).Cross(A[i] / kC_SI) ) ) / ( D * pow(1 - N.Dot(V[i] / kC_SI), 2) ) * std::exp(Exponent) * h; // FF only
        TVector3DC const ThisEw = ( ( (1 - (V[i] / kC_SI).Mag2()) * (N - (V[i] / kC_SI)) ) / ( D * D * (pow(1 - N.Dot(V[i] / kC_SI), 2)) )
            + ( N.Cross( (N - V[i] / kC_SI).Cross(A[i] / kC_SI) ) ) / ( D * pow(1 - N.Dot(V[i] / kC_SI), 2) ) ) * std::exp(Exponent) * h; // NF + FF

        SumE += ThisEw;

      }


      SumE *= C0;

      ofSpec << iw * 4.1357e-15 / k2PI << "  " <<  2 * k2PI * Current / (kh * fabs(kECharge) * kMu0 * kC_SI) *  SumE.Dot( SumE.CC() ).real()  * 1e-6  * 0.001 << std::endl;
      //ofSpec << iw * 4.1357e-15 / k2PI << "  " <<  2 * k2PI * Current / (kh * fabs(kECharge) * kMu0 * kC_SI) *  SumE.Dot( SumE.CC() ).real()  * 1e-6  * 0.001 * Obs.Mag2() << std::endl;

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

    int const NEPoints = 2000;
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

      ofSpec << iw * 4.1357e-15 / k2PI << "  " <<  8 * kPI * kPI * kEpsilon0 * kC_SI * kC_SI * Current / (kh * fabs(kECharge)) * ( SumE.Cross(SumB.CC()).Dot( TVector3DC(0, 0, 1) ) ).real() * 1e-6 * 0.001 * pow(Obs.GetZ(), 2) << std::endl;

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




  if (false) {
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
  }

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


  TBF =(TBField*) new TBField3DZRegularized(argv[1]);
  //TBF =(TBField*) new TBFieldSquareWave(0.200, 11, 0, 1.0);
  //TBF =(TBField*) new TBFieldIdeal1D(0.03, 33, 0, 1.000);
  //TBF =(TBField*) new TBFieldUniformB(0, 0, 0);



  RK4Test();

  delete TBF;

  return 0;
}

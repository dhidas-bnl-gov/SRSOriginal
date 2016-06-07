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
#include "TBFieldContainer.h"

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




//double const XStart = -1.8076921081543;
double const XStart =  -1.5;
double const XStop  =  -XStart;




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


  // Time step.  Expecting it to be constant throughout calculation
  double const DeltaT = T.GetDeltaT();


  // Number of points in the trajectory
  size_t const NTPoints = T.GetNPoints();

  // Number of points in the spectrum container
  size_t const NEPoints = S.GetNPoints();

  // Constant C0 for calculation
  double const C0 = TSRS::Qe() / (TSRS::FourPi() * TSRS::C() * TSRS::Epsilon0() * TSRS::Sqrt2Pi());

  // Constant for flux calculation at the end
  double const C2 = TSRS::FourPi() * Current / (TSRS::H() * fabs(TSRS::Qe()) * TSRS::Mu0() * TSRS::C()) * 1e-6 * 0.001;

  // Imaginary "i" and complxe 1+0i
  std::complex<double> const I(0, 1);
  std::complex<double> const One(1, 0);


  // Loop over all points in the spectrum container
  for (size_t i = 0; i != NEPoints; ++i) {

    // Angular frequency
    double const Omega = S.GetAngularFrequency(i);

    // Constant for field calculation
    std::complex<double> ICoverOmega = I * TSRS::C() / Omega;

    // Constant for calculation
    std::complex<double> const C1(0, C0 * Omega);

    // Electric field summation in frequency space
    TVector3DC SumE(0, 0, 0);

    // Loop over all points in trajectory
    for (int iT = 0; iT != NTPoints; ++iT) {

      // Particle position
      TVector3D const& X = T.GetX(iT);

      // Particle "Beta" (velocity over speed of light)
      TVector3D const& B = T.GetB(iT);

      // vector pointing from particle to observer
      TVector3D const R = ObservationPoint - X;

      // Unit vector pointing from particl to observer
      TVector3D const N = R.UnitVector();

      // Distance from particle to observer
      double const D = R.Mag();

      // Exponent for fourier transformed field
      std::complex<double> Exponent(0, Omega * (DeltaT * iT + D / TSRS::C()));

      // Sum in fourier transformed field (integral)
      SumE += (TVector3DC(B) - (N * ( One + (ICoverOmega / (D))))) / D * std::exp(Exponent);
    }

    // Multiply field by Constant C1 and time step
    SumE *= C1 * DeltaT;

    // Set the flux for this frequency / energy point
    S.SetFlux(i, C2 *  SumE.Dot( SumE.CC() ).real());
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







void CalculateFluxSurface (TParticleTrajectoryPoints const& T, TSurfacePoints const& Surface, double const Current, double const Energy_eV)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // T - Trajectory of particle
  // Surface - Surface of Observation Points
  // Current - beam current
  // Energy - beam energy in eV

  // Number of points in trajectory
  size_t const NTPoints = T.GetNPoints();

  // Time step size
  double const DeltaT = T.GetDeltaT();

  double const C0 = TSRS::Qe() / (4 * TSRS::Pi() * TSRS::C() * TSRS::Epsilon0() * TSRS::Sqrt2Pi());


  // Constant for flux calculation at the end
  double const C2 = TSRS::FourPi() * Current / (TSRS::H() * fabs(TSRS::Qe()) * TSRS::Mu0() * TSRS::C()) * 1e-6 * 0.001;

  std::complex<double> const I(0, 1);

  double const Omega = Energy_eV * TSRS::TwoPi() / 4.1357e-15;

  std::complex<double> const C1(0, C0 * Omega);

  std::ofstream of("out_flux.dat");
  of << std::scientific;

  for (size_t ip = 0; ip != Surface.GetNPoints(); ++ip) {
    if (ip % 100 == 0) {
      std::cout << ( (int) ((double) ip / (double) Surface.GetNPoints() * 100.) ) << " \% done" << std::endl;
    }

    TVector3D Obs = Surface.GetPoint(ip).GetPoint();

    TVector3DC SumE(0, 0, 0);

    for (int iT = 0; iT != NTPoints; ++iT) {
      TVector3D const& X = T.GetX(iT);
      TVector3D const& B = T.GetB(iT);
      TVector3D const& AoverC = T.GetAoverC(iT);

      TVector3D const R = Obs - X;
      TVector3D const N = R.UnitVector();
      double const D = R.Mag();
      std::complex<double> Exponent(0, -Omega * (DeltaT * iT + D / TSRS::C()));

      // TVector3DC const ThisEw = ( N.Cross( (N - B).Cross(AoverC) ) ) / ( D * pow(1 - N.Dot(B), 2) ) * std::exp(Exponent) * DeltaT; // FF only
      TVector3DC const ThisEw = ( ( (1 - (B).Mag2()) * (N - B) ) / ( D * D * (pow(1 - N.Dot(B), 2)) )
          + ( N.Cross( (N - B).Cross(AoverC) ) ) / ( D * pow(1 - N.Dot(B), 2) ) ) * std::exp(Exponent) * DeltaT; // NF + FF

      SumE += ThisEw;

    }


    SumE *= C0;

    of << Surface.GetX1(ip) << "  " << Surface.GetX2(ip) << "  " <<  C2 * SumE.Dot( SumE.CC() ).real() << std::endl;
  }

  of.close();
  exit(0);


  return;
}















void derivs(double t, double x[], double dxdt[])
{
  dxdt[0] = x[1];
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
  double x[N] = {0.0, 0.0, 0.0, 0.0, (double) XStart, (double) (BetaZ * TSRS::C())};
  double dxdt[N];


  int const NPointsForward = 20001;
  int const NPointsBack = 0; //0.5 / ((XStop - XStart) / (NPointsForward - 1));
  double const h = (XStop - XStart) / (BetaZ * TSRS::C()) / (NPointsForward - 1);

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
    x[5] = -BetaZ * TSRS::C();

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


  ParticleTrajectory.WriteToFile("del.txt");





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

    gVX.SetPoint(i, TSRS::C() * t, B.GetX());
    gVY.SetPoint(i, TSRS::C() * t, B.GetY());
    gVZ.SetPoint(i, TSRS::C() * t, B.GetZ());
    gAX.SetPoint(i, TSRS::C() * t, A.GetX());
    gAY.SetPoint(i, TSRS::C() * t, A.GetY());
    gAZ.SetPoint(i, TSRS::C() * t, A.GetZ());

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


  exit(0);





  //TVector3D Obs(0.000, 0.000, 50.7 + 0.6);
  TVector3D Obs(0.000, 0.000, 50.7 + 1.82);

  double Current = 0.500;

  TSpectrumContainer Spectrum(5000, 11500, 11900);
  CalculateSpectrumAtPoint2(ParticleTrajectory, Obs, 0.500, Spectrum);
  //CalculateSpectrumAtPoint(ParticleTrajectory, Obs, 0.500, Spectrum);
  Spectrum.SaveToFile("out_spec.dat", "");

  // Calculate Power density on a surface
  TSurfacePoints_RectangleSimple Surface("XY", 101, 101, 0.006, 0.006, 0, 0, 30, 1);
  //CalculatePowerDensitySurface(ParticleTrajectory, Surface, 0.500);

  double const Energy = 575;
  //CalculateFluxSurface(ParticleTrajectory, Surface, Current, Energy);


  // CalculateTotalPower()
  // CalculatePowerDensity()
  // CalculatePowerDensity() with gaussian convolution

  // CalculateFlux input SurfacePoint or whole Surface?
  //  How to keep data correlated...







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



  exit(0);


  // Calculate total power out
  double TotalPower = 0;
  for (int i = 0; i != NPointsForward; ++i) {
    TotalPower += ( (A[i] / TSRS::C()).Mag2() - ( (V[i] / TSRS::C()).Cross(A[i] / TSRS::C())).Mag2() ) * h;
  }
  TotalPower *= fabs(kECharge * Current) * pow(Gamma, 6) / (6 * kPI * kEpsilon0 * TSRS::C());
  std::cout << "Total Power: " << TotalPower << std::endl;

  TotalPower = 0;
  for (int i = 0; i != NPointsForward; ++i) {
    TotalPower += ( (A[i] / TSRS::C()).Mag2() - ( (V[i] / TSRS::C()).Cross(A[i] / TSRS::C())).Mag2() ) * h;
  }
  TotalPower *= fabs(kECharge * Current) * pow(Gamma, 6) / (6 * kPI * kEpsilon0 * TSRS::C());
  std::cout << "Total Power: " << TotalPower << std::endl;


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

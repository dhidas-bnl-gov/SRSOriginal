////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue May 24 11:33:19 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "SRS.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>

#include "TVector3DC.h"
#include "TBField1DZRegularized.h"
#include "TBField3DZRegularized.h"
#include "TSpectrumContainer.h"


SRS::SRS ()
{
  // Default constructor
  //fParticle.SetParticleType("");
  fCTStart = 0;
  fCTStop  = 0;
  fNPointsTrajectory = 0;
}



SRS::~SRS ()
{
  // Destructor
}



void SRS::AddMagneticField (std::string const FileName, std::string const Format, TVector3D const& Rotation, TVector3D const& Translation, std::vector<double> const& Scaling)
{
  // Add a magnetic field from a file to the field container

  std::vector<bool> HasXB(6, false);

  std::vector<int> Order(6, -1);

  std::istringstream s;
  s.str(Format);

  std::string c;
  int i = 0;
  int XDIM = 0;
  int BDIM = 0;
  while (s >> c) {


    if (c == "X") {
      Order[0] = i;
      ++XDIM;
    } else if (c == "Y") {
      Order[1] = i;
      ++XDIM;
    } else if (c == "Z") {
      Order[2] = i;
      ++XDIM;
    } else if (c == "Bx") {
      Order[3] = i;
      ++BDIM;
    } else if (c == "By") {
      Order[4] = i;
      ++BDIM;
    } else if (c == "Bz") {
      Order[5] = i;
      ++BDIM;
    } else {
      std::cerr << "ERROR: Incorrect format" << std::endl;
      throw std::invalid_argument("only excepts X Y Z Bx By Bz");
    }

    ++i;
  }


  if (XDIM > 3 || BDIM > 3) {
    std::cerr << "ERROR: spatial or B-field dimensions are too large(>3)" << std::endl;
    throw std::out_of_range("spatial or B-field dimensions are too large");
  }

  std::vector<size_t> PrefOrder;
  for (size_t j = 0; j != 6; ++j) {
    if (Order[j] != -1) {
      HasXB[i] = true;

      PrefOrder.push_back(Order[j]);
    }
  }

  // UPDATE: This needs to be fully 3D

  if (XDIM == 1) {
    if (BDIM == 1) {
      throw std::invalid_argument("Not implemented yet");
      //this->fBFieldContainer.AddField( new TBField1DZRegularized(FileName) );
      //this->fBFieldContainer.AddField( new TBField1DZRegularized(FileName, Rotation, Translation, Scaling) );
    } else if (BDIM == 2) {
      throw std::invalid_argument("Not implemented yet");
    } else if (BDIM == 3) {
      this->fBFieldContainer.AddField( new TBField3DZRegularized(FileName, Rotation, Translation, Scaling) );
    }
  } else if (XDIM == 2) {
    throw std::invalid_argument("Not implemented yet");
  } else if (XDIM == 3) {
    throw std::invalid_argument("Not implemented yet");
  } else {
    throw std::invalid_argument("Not implemented yet");
  }


  return;
}




void SRS::AddMagneticField (TBField* Field)
{
  // Add a magnetic field from a file to the field container

  this->fBFieldContainer.AddField(Field);

  return;
}




void SRS::ClearMagneticFields ()
{
  // Add a magnetic field from a file to the field container

  this->fBFieldContainer.Clear();

  return;
}




double SRS::GetBx (double const X, double const Y, double const Z) const
{
  // Return summed Bx from container
  return this->fBFieldContainer.GetBx(X, Y, Z);
}





double SRS::GetBy (double const X, double const Y, double const Z) const
{
  // Return summed By from container
  return this->fBFieldContainer.GetBy(X, Y, Z);
}





double SRS::GetBz (double const X, double const Y, double const Z) const
{
  // Return summed Bx from container
  return this->fBFieldContainer.GetBz(X, Y, Z);
}




TVector3D SRS::GetB (double const X, double const Y, double const Z) const
{
  // Return summed Bx from container
  return this->fBFieldContainer.GetB(X, Y, Z);
}




TVector3D SRS::GetB (TVector3D const& X) const
{
  // Return summed Bx from container
  return this->fBFieldContainer.GetB(X);
}




void SRS::AddParticleBeam (std::string const& Type, std::string const& Name, TVector3D const& X0, TVector3D const& V0, double const Energy_GeV, double const T0, double const Current, double const Weight, double const Charge, double const Mass)
{
  // Add a particle beam
  // Type        - The name of the particle type that you want to use
  // Name        - A user specified 'name' for this beam
  // X0          - Initial position in X,Y,Z
  // V0          - A vector pointing in the direction of the velocity of arbitrary magnitude
  // Energy_GeV  - Energy of particle beam in GeV
  // T0          - Time of initial conditions, specified in units of m (for v = c)
  // Current     - Beam current in Amperes
  // Weight      - Relative weight to give this beam when randomly sampling
  // Charge      - Charge of custom particle
  // Mass        - Mass of custom particle

  fParticleBeamContainer.AddNewParticleBeam(Type, Name, X0, V0, Energy_GeV, T0, Current, Weight, Charge, Mass);
  return;
}




TParticleBeam& SRS::GetParticleBeam (std::string const& Name)
{
  // Return a reference to the particle beam by a given name
  return fParticleBeamContainer.GetParticleBeam(Name);
}




size_t SRS::GetNParticleBeams () const
{
  // Return the number of particle beams defined
  return fParticleBeamContainer.GetNParticleBeams();
}




TParticleA SRS::GetNewParticle ()
{
  // Get a new particle.  Randomly sampled according to input beam parameters and beam weights
  return fParticleBeamContainer.GetNewParticle();
}




TParticleA const& SRS::GetCurrentParticle () const
{
  // Get a new particle.  Randomly sampled according to input beam parameters and beam weights
  return fParticle;
}




void SRS::SetNewParticle ()
{
  // Get a new particle.  Randomly sampled according to input beam parameters and beam weights.
  // Set this new particle as *the* particle in SRS fParticle
  fParticle = fParticleBeamContainer.GetNewParticle();

  return;
}




void SRS::SetNewParticle (std::string const& BeamName, std::string const& IdealOrRandom)
{
  // Get a new particle.  Randomly sampled according to input beam parameters and beam weights.
  // Set this new particle as *the* particle in SRS fParticle

  if (BeamName == "") {
    fParticle = fParticleBeamContainer.GetRandomBeam().GetNewParticle(IdealOrRandom);
  } else {
    fParticle = fParticleBeamContainer.GetParticleBeam(BeamName).GetNewParticle(IdealOrRandom);
  }

  return;
}




void SRS::ClearParticleBeams ()
{
  // Clear the contents of the particle beam container
  fParticleBeamContainer.Clear();

  return;
}




void SRS::SetNPointsTrajectory (size_t const N)
{
  // Set this number of points for any trajectory calculations
  fNPointsTrajectory = N;
  return;
}




void SRS::SetCTStart (double const X)
{
  // Set the start time in units of m (where v = c)
  fCTStart = X;
  return;
}




void SRS::SetCTStop (double const X)
{
  // Set the stop time in units of m (where v = c)
  fCTStop = X;
  return;
}




void SRS::SetCTStartStop (double const Start, double const Stop)
{
  // Set the start and stop time in units of m (where v = c)

  // If npoints is zero set this to some default value
  if (fNPointsTrajectory == 0) {
    fNPointsTrajectory = 10001 * (Stop - Start);
  }

  this->SetCTStart(Start);
  this->SetCTStop(Stop);
  return;
}




size_t SRS::GetNPointsTrajectory () const
{
  // Return the number of points being used for trajectory calculations
  return fNPointsTrajectory;
}




double SRS::GetCTStart () const
{
  // Return the start time in units of m (where v = c)
  return fCTStart;
}




double SRS::GetCTStop () const
{
  // Return the stop time in units of m (where v = c)
  return fCTStop;
}




void SRS::CalculateTrajectory ()
{
  // Function to calculate the particle trajectory of the member particle fParticle

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    this->SetNewParticle();
  }

  this->CalculateTrajectory(fParticle);

  return;
}





void SRS::CalculateTrajectory (TParticleA& P)
{
  // Function to calculate the particle trajectory given initial conditions.
  // This function uses the internal Trajectory member to store results

  // UPDATE: Check this hits correct points in trajectory and correct number of points
  // UPDATE: Do a detailed check of forward and back prop of X, V, and A

  // Check that CTStart is not after T0 of particle
  if (this->GetCTStart() > P.GetT0()) {
    std::cerr << "ERROR: start time is greater than T0" << std::endl;
    throw std::out_of_range("start time is greater than T0");
  }

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (P.GetType() == "") {
    throw std::out_of_range("particle not initialized.  make sure you have a particle or beam defined");
  }

  // Clear any current trajectory
  P.GetTrajectory().Clear();


  // Calculate the total DeltaT in seconds
  double const DeltaT = ((this->GetCTStop() - this->GetCTStart()) / TSRS::C() / (fNPointsTrajectory - 1));


  // The number of points in the forward and backward direction
  // UPDATE: Check the details of these numbers
  size_t NPointsForward  = 1 + (this->GetCTStop() - P.GetT0()) / TSRS::C() / DeltaT;
  size_t NPointsBackward = (P.GetT0() - this->GetCTStart()) / TSRS::C() / DeltaT;


  // Number of dimensions of the array to be sent for RK calculation
  int const N = 6;

  // Arrays to be sent for RK calculation
  double x[N];
  double dxdt[N];


  // Initial conditions for the forward propogation
  x[0] = P.GetX0().GetX();
  x[1] = P.GetB0().GetX() * TSRS::C();
  x[2] = P.GetX0().GetY();
  x[3] = P.GetB0().GetY() * TSRS::C();
  x[4] = P.GetX0().GetZ();
  x[5] = P.GetB0().GetZ() * TSRS::C();


  // Grap the particle trajectory object
  TParticleTrajectoryPoints& ParticleTrajectory = P.GetTrajectory();

  // Set delta T for the trajectory
  ParticleTrajectory.SetDeltaT(DeltaT);

  // Loop over points in the forward direction
  for (int i = 0; i != NPointsForward; ++i) {

    // This time
    double t = P.GetT0() + DeltaT * i;


    // Add this point to the trajectory
    ParticleTrajectory.AddPoint(x[0], x[2], x[4], x[1] / TSRS::C(), x[3] / TSRS::C(), x[5] / TSRS::C(), dxdt[1] / TSRS::C(), dxdt[3] / TSRS::C(), dxdt[5] / TSRS::C());

    // Propogate
    this->Derivatives(t, x, dxdt, P);
    RK4(x, dxdt, N, t, DeltaT, x, P);
  }

  // Reverse trajectory elements for backward propogation
  ParticleTrajectory.ReverseArrays();

  // Set initial conditions for propogating backwards
  x[0] =  P.GetX0().GetX();
  x[1] =  P.GetB0().GetX() * TSRS::C();
  x[2] =  P.GetX0().GetY();
  x[3] =  P.GetB0().GetY() * TSRS::C();
  x[4] =  P.GetX0().GetZ();
  x[5] =  P.GetB0().GetZ() * TSRS::C();

  // Reverse time
  double const DeltaTReversed = -DeltaT;

  // Loop over all points "before" the initial point
  for (int i = 0; i != NPointsBackward; ++i) {

    // This time
    double t = P.GetT0() + DeltaTReversed * (i + 1);

    // Propogate backward in time!
    this->Derivatives(t, x, dxdt, P);
    RK4(x, dxdt, N, t, DeltaTReversed, x, P);

    // Add the point to the trajectory
    ParticleTrajectory.AddPoint(x[0], x[2], x[4], x[1] / TSRS::C(), x[3] / TSRS::C(), x[5] / TSRS::C(), dxdt[1] / TSRS::C(), dxdt[3] / TSRS::C(), dxdt[5] / TSRS::C());
  }

  // Re-Reverse the trajectory to be in the proper time order
  ParticleTrajectory.ReverseArrays();

  return;
}




TParticleTrajectoryPoints const& SRS::GetTrajectory ()
{
  // Get the trajectory for *the* current particle in fParticle

  return fParticle.GetTrajectory();
}











void SRS::Derivatives (double t, double x[], double dxdt[], TParticleA const& P)
{
  // This is a second order differential equation.  It does not account for the loss in energy due to
  // radiation.  Although 't' is not used it would be easy to implement a time dependent field

  // The values correspond to:
  // x[0] - x
  // x[1] - Vx
  // x[2] - y
  // x[3] - Vy
  // x[4] - z
  // x[5] - Vz

  // BField at this point
  TVector3D const B = this->GetB(x[0], x[2], x[4]);

  dxdt[0] = x[1];
  dxdt[1] = P.GetQoverMGamma() * (-x[5] * B.GetY() + x[3] * B.GetZ());
  dxdt[2] = x[3];                                                                                            
  dxdt[3] = P.GetQoverMGamma() * ( x[5] * B.GetX() - x[1] * B.GetZ());
  dxdt[4] = x[5];                                                                                            
  dxdt[5] = P.GetQoverMGamma() * ( x[1] * B.GetY() - x[3] * B.GetX());

  return;
}



void SRS::RK4 (double y[], double dydx[], int n, double x, double h, double yout[], TParticleA const& P)
{
  // Runge-Kutta 4th order method propogation
  // UPDATE: syntax

  int i;
  double xh, hh, h6;
  double *dym, *dyt, *yt;

  dym = new double[n];
  dyt = new double[n];
  yt  = new double[n];

  hh = h * 0.5;
  h6 = h / 6.0;
  xh = x + hh;

  for (i = 0; i < n; ++i) {
    yt[i] = y[i] + hh * dydx[i];
  }

  this->Derivatives(xh, yt, dyt, P);

  for (i = 0; i < n; ++i) {
    yt[i] = y[i] + hh * dyt[i];
  }

  this->Derivatives(xh, yt, dym, P);

  for (i = 0; i < n; ++i) {
    yt[i] = y[i] + h * dym[i];
    dym[i] += dyt[i];
  }

  this->Derivatives(x + h, yt, dyt, P);

  for (i = 0; i < n; ++i) {
    yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
  }
  
  delete [] dym;
  delete [] dyt;
  delete [] yt;

  return;
}




void SRS::CalculateSpectrum (TVector3D const& ObservationPoint, TSpectrumContainer& Spectrum, double const Weight)
{
  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  this->CalculateSpectrum(fParticle, ObservationPoint, Spectrum, Weight);
  return;
}




void SRS::CalculateSpectrum (TParticleA& Particle, TVector3D const& ObservationPoint, TSpectrumContainer& Spectrum, double const Weight)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  // Save this in the spectrum container.
  //
  // Particle - the Particle.. with a Trajectory structure hopefully
  // ObservationPoint - Observation Point
  // Spectrum - Spectrum container

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    throw std::out_of_range("no particle defined");
  }

  // Calculate trajectory
  this->CalculateTrajectory(Particle);

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();


  // Time step.  Expecting it to be constant throughout calculation
  double const DeltaT = T.GetDeltaT();


  // Number of points in the trajectory
  size_t const NTPoints = T.GetNPoints();

  if (NTPoints < 1) {
    throw std::length_error("no points in trajectory.  Is particle or beam defined?");
  }

  // Number of points in the spectrum container
  size_t const NEPoints = Spectrum.GetNPoints();

  // Constant C0 for calculation
  double const C0 = Particle.GetQ() / (TSRS::FourPi() * TSRS::C() * TSRS::Epsilon0() * TSRS::Sqrt2Pi());

  // Constant for flux calculation at the end
  double const C2 = TSRS::FourPi() * Particle.GetCurrent() / (TSRS::H() * fabs(Particle.GetQ()) * TSRS::Mu0() * TSRS::C()) * 1e-6 * 0.001;

  // Imaginary "i" and complxe 1+0i
  std::complex<double> const I(0, 1);
  std::complex<double> const One(1, 0);


  // Loop over all points in the spectrum container
  for (size_t i = 0; i != NEPoints; ++i) {

    // Angular frequency
    double const Omega = Spectrum.GetAngularFrequency(i);

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

      // Vector pointing from particle to observer
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
    Spectrum.AddToFlux(i, C2 *  SumE.Dot( SumE.CC() ).real() * Weight);
  }



  return;
}




void SRS::CalculateSpectrum (TParticleA& Particle, TVector3D const& ObservationPoint, double const EStart, double const EStop, size_t const N, std::string const& OutFilename)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // Particle - the Particle.. with a Trajectory structure hopefully
  // ObservationPoint - Observation Point

  // A spectrum object for return values
  fSpectrum.Init(N, EStart, EStop);

  this->CalculateSpectrum(Particle, ObservationPoint, fSpectrum);

  return;
}




void SRS::CalculateSpectrum (TVector3D const& ObservationPoint, double const EStart, double const EStop, size_t const N)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // ObservationPoint - Observation Point

  // A spectrum object for return values
  fSpectrum.Init(N, EStart, EStop);

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  this->CalculateSpectrum(fParticle, ObservationPoint, fSpectrum);

  return;
}




void SRS::CalculateSpectrum (TVector3D const& ObservationPoint, std::vector<double> const& V)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // Particle - the Particle.. with a Trajectory structure hopefully
  // ObservationPoint - Observation Point

  // A spectrum object for return values
  fSpectrum.Init(V);

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  this->CalculateSpectrum(fParticle, ObservationPoint, fSpectrum);

  return;
}




TSpectrumContainer const& SRS::GetSpectrum () const
{
  return fSpectrum;
}




void SRS::CalculatePowerDensity (TParticleA& Particle, TSurfacePoints const& Surface, int const Dimension, bool const Directional, std::string const& OutFileName)
{
  T3DScalarContainer PowerDensityContainer;
  this->CalculatePowerDensity(Particle, Surface, PowerDensityContainer, Dimension, Directional, OutFileName);

  return;
}




void SRS::CalculatePowerDensity (TParticleA& Particle, TSurfacePoints const& Surface, T3DScalarContainer& PowerDensityContainer, int const Dimension, bool const Directional, std::string const& OutFileName)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // Particle - Particle, contains trajectory (or if not, calculate it)
  // Surface - Observation Point


  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    throw std::out_of_range("no particle defined");
  }

  // Are we writing to a file?
  bool const WriteToFile = OutFileName != "" ? true : false;

  // Calculate the trajectory from scratch
  this->CalculateTrajectory(Particle);

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();

  // Number of points in Trajectory
  size_t const NTPoints = T.GetNPoints();

  // Timestep from trajectory
  double const DeltaT = T.GetDeltaT();

  // Variables for parts of the numerator and denominator of power density equation
  TVector3D Numerator;
  double Denominator;

  // If writing to a file, open it and set to scientific output
  std::ofstream of;
  if (WriteToFile) {
    of.open(OutFileName.c_str());
    if (!of.is_open()) {
      throw std::ofstream::failure("cannot open output file");
    }
    of << std::scientific;
  }
  //std::cout << "Directional: " << Directional << std::endl;

  // Loop over all points in the given surface
  for (size_t io = 0; io < Surface.GetNPoints(); ++io) {

    // Get the observation point (on the surface, and its "normal"
    TVector3D Obs = Surface.GetPoint(io).GetPoint();
    TVector3D Normal = Surface.GetPoint(io).GetNormal();


    // For summing power contributions
    double Sum = 0;

    // Loop over all points in the trajectory
    for (int iT = 0; iT != NTPoints ; ++iT) {

      // Get current position, Beta, and Acceleration(over c)
      TVector3D const& X = T.GetX(iT);
      TVector3D const& B = T.GetB(iT);
      TVector3D const& AoverC = T.GetAoverC(iT);

      // Define the three normal vectors.  N1 is in the direction of propogation,
      // N2 and N3 are in a plane perpendicular to N1
      TVector3D const N1 = (Obs - X).UnitVector();
      TVector3D const N2 = N1.Orthogonal().UnitVector();
      TVector3D const N3 = N1.Cross(N2).UnitVector();

      // For computing non-normally incidence
      double const N1DotNormal = N1.Dot(Normal);
      // UPDATE: URGENT: Check this
      if (Directional && N1DotNormal <= 0) {
        continue;
      }

      // Compute Numerator and denominator
      Numerator = N1.Cross( ( (N1 - B).Cross((AoverC)) ) );
      Denominator = pow(1 - (B).Dot(N1), 5);

      // Add contributions from both N2 and N3
      Sum += pow(Numerator.Dot(N2), 2) / Denominator / (Obs - X).Mag2() * N1DotNormal;
      Sum += pow(Numerator.Dot(N3), 2) / Denominator / (Obs - X).Mag2() * N1DotNormal;
    }

    // Undulators, Wigglers and their applications, p42
    Sum *= fabs(Particle.GetQ() * Particle.GetCurrent()) / (16 * TSRS::Pi2() * TSRS::Epsilon0() * TSRS::C()) * DeltaT;

    Sum /= 1e6; // m^2 to mm^2

    // If you don't care about the direction of the normal vector
    // UPDATE: Check
    if (!Directional) {
      if (Sum < 0) {
        Sum *= -1;
      }
    }

    if (Dimension == 2) {
      if (WriteToFile) {
        of << Surface.GetX1(io) << " " << Surface.GetX2(io) << " " << Sum << "\n";
      } else {
        PowerDensityContainer.AddPoint(TVector3D(Surface.GetX1(io), Surface.GetX2(io), 0), Sum);
      }
    } else if (Dimension == 3) {
      if (WriteToFile) {
        of << Obs.GetX() << " " << Obs.GetY() << " " << Obs.GetZ() << " " << Sum << "\n";
      } else {
        PowerDensityContainer.AddPoint(Obs, Sum);
      }
    } else {
      throw std::out_of_range("incorrect dimensions");
    }

  }

  if (WriteToFile) {
    of.close();
  }

  return;
}




void SRS::CalculatePowerDensity (TSurfacePoints const& Surface, T3DScalarContainer& PowerDensityContainer, int const Dimension, bool const Directional, std::string const& OutFileName)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // Surface - Observation Point

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  this->CalculatePowerDensity(fParticle, Surface, PowerDensityContainer, Dimension, Directional, OutFileName);

  return;
}





double SRS::CalculateTotalPower ()
{
  // UPDATE: commet

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  return this->CalculateTotalPower(fParticle);
}





double SRS::CalculateTotalPower (TParticleA& Particle)
{
  // Calculate total power out

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    throw std::out_of_range("no particle defined");
  }

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();

  // Calculate Trajectory from scratch
  this->CalculateTrajectory(Particle);

  // Number of points in Trajectory
  size_t const NTPoints = T.GetNPoints();

  // Timestep from trajectory
  double const DeltaT = T.GetDeltaT();

  // For summing total power
  double TotalPower = 0;


  // Loop over all points in trajectory
  for (int i = 0; i != NTPoints; ++i) {
    TVector3D const& B = T.GetB(i);
    TVector3D const& AoverC = T.GetAoverC(i);

    TotalPower += (AoverC.Mag2() - (B.Cross(AoverC)).Mag2()) * DeltaT;

  }

  // Undulators, Wigglers and their applications, p42
  TotalPower *= fabs(Particle.GetQ() * Particle.GetCurrent()) * pow(Particle.GetGamma(), 6) / (6 * TSRS::Pi() * TSRS::Epsilon0() * TSRS::C());


  return TotalPower;
}





void SRS::CalculateFlux (TParticleA& Particle, TSurfacePoints const& Surface, double const Energy_eV, std::string const& OutFileName)
{
  T3DScalarContainer FluxContainer;
  this->CalculateFlux(Particle, Surface, Energy_eV, FluxContainer);

  return;
}








void SRS::CalculateFlux (TParticleA& Particle, TSurfacePoints const& Surface, double const Energy_eV, T3DScalarContainer& FluxContainer, std::string const& OutFileName)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // T - Trajectory of particle
  // Surface - Surface of Observation Points
  // Current - beam current
  // Energy - beam energy in eV

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    throw std::out_of_range("no particle defined");
  }

  int const Dimension = 2;

  // Write to a file?
  bool const WriteToFile = OutFileName != "" ? true : false;


  // Calculate trajectory
  this->CalculateTrajectory(Particle);

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();

  // Number of points in trajectory
  size_t const NTPoints = T.GetNPoints();

  // Time step size
  double const DeltaT = T.GetDeltaT();

  // Constand C0
  double const C0 = Particle.GetQ() / (4 * TSRS::Pi() * TSRS::C() * TSRS::Epsilon0() * TSRS::Sqrt2Pi());


  // Constant for flux calculation at the end
  double const C2 = TSRS::FourPi() * Particle.GetCurrent() / (TSRS::H() * fabs(Particle.GetQ()) * TSRS::Mu0() * TSRS::C()) * 1e-6 * 0.001;

  // Complex number i
  std::complex<double> const I(0, 1);

  // UPDATE: use common function
  double const Omega = Energy_eV * TSRS::TwoPi() / 4.1357e-15;

  // Constant C1 (complex)
  //std::complex<double> const C1(0, C0 * Omega);

  // If writing to a file, open it and set to scientific output
  std::ofstream of;
  if (WriteToFile) {
    of.open(OutFileName.c_str());
    if (!of.is_open()) {
      throw std::ofstream::failure("cannot open output file");
    }
    of << std::scientific;
  }

  // Loop over all surface points
  for (size_t ip = 0; ip != Surface.GetNPoints(); ++ip) {

    // Observation point
    TVector3D Obs = Surface.GetPoint(ip).GetPoint();

    // Sum E-field
    TVector3DC SumE(0, 0, 0);

    // Loop over trajectory points
    for (int iT = 0; iT != NTPoints; ++iT) {

      // Get position, Beta, and Acceleration (over c)
      TVector3D const& X = T.GetX(iT);
      TVector3D const& B = T.GetB(iT);
      TVector3D const& AoverC = T.GetAoverC(iT);

      // Define R and unit vector in direction of R, and D (distance to observer)
      TVector3D const R = Obs - X;
      TVector3D const N = R.UnitVector();
      double const D = R.Mag();

      // Exponent in transformed field
      std::complex<double> Exponent(0, -Omega * (DeltaT * iT + D / TSRS::C()));

      // UPDATE: with better field avoiding A
      // TVector3DC const ThisEw = ( N.Cross( (N - B).Cross(AoverC) ) ) / ( D * pow(1 - N.Dot(B), 2) ) * std::exp(Exponent) * DeltaT; // FF only

      TVector3DC const ThisEw = ( ( (1 - (B).Mag2()) * (N - B) ) / ( D * D * (pow(1 - N.Dot(B), 2)) )
          + ( N.Cross( (N - B).Cross(AoverC) ) ) / ( D * pow(1 - N.Dot(B), 2) ) ) * std::exp(Exponent) * DeltaT; // NF + FF

      // Add this contribution
      SumE += ThisEw;

    }

    // Multiply by constant factor
    SumE *= C0;

    //of << Surface.GetX1(ip) << "  " << Surface.GetX2(ip) << "  " <<  C2 * SumE.Dot( SumE.CC() ).real() << std::endl;
    //FluxContainer.AddPoint(Obs, C2 * SumE.Dot( SumE.CC() ).real());

    double const ThisFlux = C2 * SumE.Dot( SumE.CC() ).real();


    // If you don't care about the direction of the normal vector
    //if (!Directional) {
    //  if (Sum < 0) {
    //    Sum *= -1;
    //  }
    //}

    if (Dimension == 2) {
      if (WriteToFile) {
        of << Surface.GetX1(ip) << " " << Surface.GetX2(ip) << " " << ThisFlux << "\n";
      } else {
        FluxContainer.AddPoint(TVector3D(Surface.GetX1(ip), Surface.GetX2(ip), 0), ThisFlux);
      }
    } else if (Dimension == 3) {
      if (WriteToFile) {
        of << Obs.GetX() << " " << Obs.GetY() << " " << Obs.GetZ() << " " << ThisFlux << "\n";
      } else {
        FluxContainer.AddPoint(Obs, ThisFlux);
      }
    } else {
      throw std::out_of_range("incorrect dimensions");
    }

  }

  if (WriteToFile) {
    of.close();
  }

  return;
}








void SRS::CalculateFlux (TSurfacePoints const& Surface, double const Energy_eV, T3DScalarContainer& FluxContainer, std::string const& OutFileName)
{
  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  this->CalculateFlux(fParticle, Surface, Energy_eV, FluxContainer);
  return;
}





void SRS::CalculateElectricFieldTimeDomain (TVector3D const& Observer, T3DScalarContainer& XYZT)
{
  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  this->CalculateElectricFieldTimeDomain(Observer, XYZT, fParticle);
  return;
}





void SRS::CalculateElectricFieldTimeDomain (TVector3D const& Observer, T3DScalarContainer& XYZT,  TParticleA& Particle)
{
  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    throw std::out_of_range("no beam defined");
  }

  // Calculate Trajectory
  this->CalculateTrajectory(Particle);

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();
  size_t const NTPoints = T.GetNPoints();

  double const DeltaT = T.GetDeltaT();

  double const C0 = Particle.GetQ() / (TSRS::FourPi() * TSRS::Epsilon0());
  // Loop over trajectory points
  for (size_t iT = 0; iT != NTPoints; ++iT) {

    // Get position, Beta, and Acceleration (over c)
    TVector3D const& X = T.GetX(iT);
    TVector3D const& B = T.GetB(iT);
    TVector3D const& AoverC = T.GetAoverC(iT);

    // Define R and unit vector in direction of R, and D (distance to observer)
    TVector3D const R = Observer - X;
    TVector3D const N = R.UnitVector();
    double const D = R.Mag();

    // Synchrotron Radiation, Philip john Duke, p74

    double    const      Mult = C0 * pow(1.0 / (1.0 - N.Dot(B)), 3);
    TVector3D const NearField = ((1.0 - B.Mag2() ) * (N - B)) / R.Mag2();
    TVector3D const  FarField = (1.0 / TSRS::C()) * (N.Cross(  (N - B).Cross(AoverC))  ) / R.Mag();

    TVector3D const EField = Mult * (NearField + FarField);
    double    const Time = iT * DeltaT + D / TSRS::C();

    XYZT.AddPoint(EField, Time);
  }

  //return TVector3D(Mult * (NearField + FarField));
  return;
}

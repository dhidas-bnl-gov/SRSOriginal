////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Jun  8 17:18:39 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TParticleBeam.h"

#include "TSRS.h"

#include <cmath>


TParticleBeam::TParticleBeam ()
{
  // Default constructor
}




TParticleBeam::TParticleBeam (std::string const& ParticleType)
{
  // Constructor given a particle type.
  // Sets the current to the single charge / s

  this->SetParticleType(ParticleType);

  this->SetCurrent(this->GetQ());
}




TParticleBeam::TParticleBeam (std::string const& ParticleType, double const Energy, double const Current)
{
  // Constructor given a particle type.

  this->SetParticleType(ParticleType);
  fE0 = Energy;

  this->SetCurrent(Current);
}




TParticleBeam::TParticleBeam (std::string const& ParticleType, TVector3D const& X0, TVector3D const& D0, double const Energy, double const Current)
{
  // Constructor given a particle type.
  // Sets the initial time to 0

  this->SetParticleType(ParticleType);

  fX0 = X0;
  fU0 = D0.UnitVector();
  fE0 = Energy;
  fT0 = 0;

  this->SetCurrent(Current);
}




TParticleBeam::TParticleBeam (std::string const& ParticleType, TVector3D const& X0, TVector3D const& D0, double const Energy, double const T0, double const Current, double const Charge, double const Mass)
{
  std::cout << "current " << Current << std::endl;
  // Constructor given a particle type.

  if (ParticleType == "custom") {
    this->SetParticleTypeCustom(ParticleType, Charge, Mass);
  } else {
    this->SetParticleType(ParticleType);
  }

  fX0 = X0;
  fU0 = D0.UnitVector();
  fE0 = Energy;
  fT0 = T0;

  this->SetCurrent(Current);
}




TParticleBeam::TParticleBeam (std::string const& ParticleType, TVector3D const& X0, TVector3D const& D0, double const Energy, double const T0, double const Current, TTwiss const& Twiss)
{
  // Constructor given a particle type.

  this->SetParticleType(ParticleType);

  fX0 = X0;
  fU0 = D0.UnitVector();
  fE0 = Energy;
  fT0 = T0;

  this->SetCurrent(Current);

  fTwiss = Twiss;
}




TParticleBeam::~TParticleBeam ()
{
  // Destruction!
}




void TParticleBeam::SetInitialConditions (double const X, double const Y, double const Z, double const Dx, double const Dy, double const Dz, double const E0, double const T)
{
  // This function takes the initial position, initial direction (normalized or not), initial time, and Energy of a particle

  this->SetInitialConditions( TVector3D(X, Y, Z), TVector3D(Dx, Dy, Dz), E0, T );
  return;
}




void TParticleBeam::SetInitialConditions (TVector3D const& X, TVector3D const& D, double const E0, double const T0)
{
  // Set the initial conditions variables for this particle
  // X - Initial position
  // D - 3D vector of the initial direction of velocity (Magnitude is arbitrary)
  // E0 - The initial energy
  // T0 - The initial time (in [m])

  // Rescale the direction to a unit vector

  this->fX0 = X;
  this->fU0 = D.UnitVector();
  this->fE0 = E0;
  this->fT0 = T0;

  return;
}




TVector3D const& TParticleBeam::GetX0 () const
{
  // Return const reference to the initial position
  return fX0;
}




TVector3D const& TParticleBeam::GetU0 () const
{
  // Get unit vector in the direction of initial velocity
  return fU0;
}




double TParticleBeam::GetE0 () const
{
  // Return initial energy
  return fE0;
}




double TParticleBeam::GetT0 () const
{
  // Return initial time
  return fT0;
}




TParticleA TParticleBeam::GetNewParticle ()
{
  // Intended to get you a new random particle based on the beam parameters
  // given.

  // UPDATE: not just below, but all
  this->GetTrajectory().Clear();

  // UPDATE: Needs rand for twiss, or other beam configurations...
  // UPDATE: Could also take a python function


  double const Gamma = this->GetE0() / TSRS::kgToGeV(this->GetM());
  double const Beta = sqrt(1.0 - 1.0 / (Gamma * Gamma));

  TVector3D XNew;
  TVector3D BNew;
  double    TNew;
  double    ENew; // correlated with BNew, not sure how to handle this yet


  TVector3D Beta0 = this->GetU0() * Beta;

  TParticleA NewParticle = (TParticleA) *this;
  NewParticle.SetInitialParticleConditions(this->GetX0(), Beta0, this->GetT0());


  return NewParticle;
}

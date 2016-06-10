////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Jun  8 14:35:28 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TParticleA.h"
#include "TSRS.h"

#include <algorithm>
#include <cmath>



TParticleA::TParticleA ()
{
}




TParticleA::TParticleA (std::string const& Type)
{
  this->SetParticleType(Type);
}




TParticleA::TParticleA (std::string const& Type, TVector3D const& X0, TVector3D const& B0, double const T0)
{
  this->SetParticleType(Type);
  this->SetX0(X0);
  this->SetB0(B0);
  this->SetT0(T0);

  SetGamma();
}




TParticleA::~TParticleA ()
{
}




void TParticleA::SetParticleType (std::string const& Type)
{
  std::string type = Type;
  std::transform(type.begin(), type.end(), type.begin(), ::tolower);

  // So some stuff

  // Leptons first
  if (type == "electron" || type == "anti-positron") {
    this->SetQM( -TSRS::Qe(), TSRS::Me() );
  } else if (type == "positron" || type == "anti-electron") {
    this->SetQM(  TSRS::Qe(), TSRS::Me() );
  } else if (type == "muon") {
    this->SetQM( -TSRS::Qe(), TSRS::Me() );
  } else if (type == "anti-muon") {
    this->SetQM(  TSRS::Qe(), TSRS::Me() );


  } else if (type == "proton") {
    this->SetQM(  TSRS::Qe(), TSRS::Me() );
  } else if (type == "anti-proton") {
    this->SetQM( -TSRS::Qe(), TSRS::Me() );
  } else if (type == "neutron") {
    this->SetQM(           0, TSRS::Me() );
  } else {
    throw;
  }

  return;
}




void TParticleA::SetParticleTypeFromPDGID (int const ID)
{
  // So some stuff
  throw;
  return;
}




void TParticleA::SetQ (double const Q)
{
  fQ = Q;

  SetQoverMGamma();

  return;
}




void TParticleA::SetM (double const M)
{
  fM = M;

  SetQoverMGamma();

  return;
}




void TParticleA::SetQM (double const Q, double const M)
{
  fQ = Q;
  fM = M;

  SetQoverMGamma();

  return;
}




void TParticleA::SetX0 (TVector3D const& X0)
{
  fX0 = X0;
  return;
}




void TParticleA::SetB0 (TVector3D const& B0)
{
  fB0 = B0;

  SetGamma();

  return;
}




void TParticleA::SetT0 (double const T0)
{
  fT0 = T0;
  return;
}




double TParticleA::GetQ () const
{
  return fQ;
}




double TParticleA::GetM () const
{
  return fM;
}




double TParticleA::GetGamma () const
{
  return fGamma;
}




double TParticleA::GetQoverMGamma () const
{
  return fQoverMGamma;
}




TVector3D const& TParticleA::GetX0 () const
{
  return fX0;
}




TVector3D const& TParticleA::GetB0 () const
{
  return fB0;
}




double TParticleA::GetT0 () const
{
  return fT0;
}




void TParticleA::SetInitialParticleConditions (TVector3D const& X0, TVector3D const& B0, double const T0)
{
  // Set the initial conditions for this particle
  fX0 = X0;
  fB0 = B0;
  fT0 = T0;

  SetGamma();

  return;
}




TParticleTrajectoryPoints& TParticleA::GetTrajectory ()
{
  return fTrajectory;
}




void TParticleA::SetGamma ()
{

  double const Beta2 = fB0.Mag2();
  if (Beta2 == 1) {
    return;
  }

  fGamma = 1.0 / sqrt(1.0 - fB0.Mag2());

  SetQoverMGamma();


  return;
}




void TParticleA::SetQoverMGamma ()
{
  if (GetM() == 0 || GetGamma() == 0) {
    return;
  }

  fQoverMGamma = GetQ() / GetM() / GetGamma();

  return;
}

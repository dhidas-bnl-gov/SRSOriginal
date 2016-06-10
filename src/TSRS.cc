////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue May 24 11:33:19 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TSRS.h"

#include "TBField3DZRegularized.h"


SRS::SRS ()
{
  // Default constructor
}



SRS::~SRS ()
{
  // Destructor
}



void SRS::AddMagneticField (std::string const FileName, std::string const Format, double const X0, double const Y0, double const Z0)
{
  // Add a magnetic field from a file to the field container

  this->fBFieldContainer.AddField( new TBField3DZRegularized(FileName) );

  return;
}




double SRS::GetBx (double const X, double const Y, double const Z) const
{
  return this->fBFieldContainer.GetBx(X, Y, Z);
}





double SRS::GetBy (double const X, double const Y, double const Z) const
{
  return this->fBFieldContainer.GetBy(X, Y, Z);
}





double SRS::GetBz (double const X, double const Y, double const Z) const
{
  return this->fBFieldContainer.GetBz(X, Y, Z);
}




void SRS::AddParticleBeam (std::string const& Type, std::string const& Name, double const X0, double const Y0, double const Z0, double const DX0, double const DY0, double const DZ0, double const Energy, double const T0, double const Current, double const Weight)
{
  fParticleBeamContainer.AddNewParticleBeam(Type, Name, TVector3D(X0, Y0, Z0), TVector3D(DX0, DY0, DZ0), Energy, T0, Current, Weight);
  return;
}




TParticleBeam& SRS::GetParticleBeam (std::string const& Name)
{
  return fParticleBeamContainer.GetParticleBeam(Name);
}




TParticleA SRS::GetNewParticle ()
{
  return fParticleBeamContainer.GetNewParticle();
}

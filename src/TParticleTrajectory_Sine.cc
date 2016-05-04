////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Apr  1 08:23:25 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TParticleTrajectory_Sine.h"


#include <cmath>


TParticleTrajectory_Sine::TParticleTrajectory_Sine (double const& Wavelength, double const& Amplitude, double const& Phase)
{
  // Constructor
  fWavelength = Wavelength;
  fAmplitude = Amplitude;
  fPhase = Phase;
}




TParticleTrajectory_Sine::~TParticleTrajectory_Sine ()
{
  // Destructor!
}




TVector3D TParticleTrajectory_Sine::GetPosition (double const& Time) const
{
  return TVector3D( fAmplitude * sin(2 * 3.14 * Time / fWavelength + fPhase), 3e8 * Time, 0);
}





TVector3D TParticleTrajectory_Sine::GetVelocity (double const& Time) const
{
  return TVector3D( (2 * 3.14 / fWavelength) * fAmplitude * cos(2 * 3.14 * Time / fWavelength + fPhase), 3e8, 0);
}




TVector3D TParticleTrajectory_Sine::GetAcceleration (double const& Time) const
{
  return TVector3D( -pow(2 * 3.14 / fWavelength, 2) * fAmplitude * sin(2 * 3.14 * Time / fWavelength + fPhase), 0, 0);
}





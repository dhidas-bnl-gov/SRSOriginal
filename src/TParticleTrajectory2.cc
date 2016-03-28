////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Mar 25 17:18:46 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TParticleTrajectory2.h"



TParticleTrajectory2::TParticleTrajectory2 ()
{
  // Default constructor
  // If you use this you need to set the following variables by hand.. I advise against it.

  fDeltaTime = 0;
  fStartTime = 0;
}




TParticleTrajectory2::TParticleTrajectory2 (double const& StartTime, double const& DeltaTime)
{
  // This is the contstrctur you should use.  Without these variables set the calculations
  // won't work

  fStartTime = StartTime;
  fDeltaTime = DeltaTime;
}




TParticleTrajectory2::~TParticleTrajectory2 ()
{
  // Destructor
  // This object owns the memory of fX, fY, and fZ and they must be deleted.

  for (std::vector<TVector3D*>::iterator it = fX.begin(); it != fX.end(); ++it) {
    delete *it;
  }
  for (std::vector<TVector3D*>::iterator it = fV.begin(); it != fV.end(); ++it) {
    delete *it;
  }
  for (std::vector<TVector3D*>::iterator it = fA.begin(); it != fA.end(); ++it) {
    delete *it;
  }
}




void TParticleTrajectory2::AddPoint(double const& X, double const& Y, double const& Z, double const& Vx, double const& Vy, double const& Vz, double const& Ax, double const& Ay, double const& Az)
{
  // Add a point in the forward time direction.  This assumes the fDeltaTime time step

  fX.push_back(new TVector3D(X, Y, Z));
  fV.push_back(new TVector3D(Vx, Vy, Vz));
  fA.push_back(new TVector3D(Ax, Ay, Az));

  return;
}




TVector3D TParticleTrajectory2::GetXAtTime (double const& Time) const
{
  // Get the position TVector3D at any time (within the allowed time range)

  size_t i = GetFirstBinForTime(Time);
  return TVector3D( *fX[i] + (*fX[i+1] - *fX[i]) * (Time - (fStartTime + fDeltaTime * (float) i)));
}




TVector3D TParticleTrajectory2::GetVAtTime (double const& Time) const
{
  // Get the velocity TVector3D at any time (within the allowed time range)

  size_t i = GetFirstBinForTime(Time);
  return TVector3D( *fV[i] + (*fV[i+1] - *fV[i]) * (Time - (fStartTime + fDeltaTime * (float) i)));
}




TVector3D TParticleTrajectory2::GetAAtTime (double const& Time) const
{
  // Get the acceleration TVector3D at any time (within the allowed time range)

  size_t i = GetFirstBinForTime(Time);
  return TVector3D( *fA[i] + (*fA[i+1] - *fA[i]) * (Time - (fStartTime + fDeltaTime * (float) i)));
}




size_t TParticleTrajectory2::GetFirstBinForTime (double const& Time) const
{
  // Get the time bin index just before the given time

  if (Time < fStartTime) {
    throw;
  }

  if (Time > fStartTime + fDeltaTime * ((double) (fX.size() - 1))) {
    throw;
  }

  return (size_t) ((Time - fStartTime) / fDeltaTime);
}

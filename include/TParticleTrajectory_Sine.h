#ifndef GUARD_TParticleTrajectory_Sine_h
#define GUARD_TParticleTrajectory_Sine_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Apr  1 08:23:25 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TParticleTrajectory.h"


class TParticleTrajectory_Sine : public TParticleTrajectory
{
  public:
    TParticleTrajectory_Sine (double const&, double const&, double const&);
    ~TParticleTrajectory_Sine ();
    TVector3D GetPosition (double const&) const;
    TVector3D GetVelocity (double const&) const;
    TVector3D GetAcceleration (double const&) const;

  private:
    double fWavelength;
    double fAmplitude;
    double fPhase;

};













#endif


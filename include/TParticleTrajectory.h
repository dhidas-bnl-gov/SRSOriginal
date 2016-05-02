#ifndef GUARD_TParticleTrajectory_h
#define GUARD_TParticleTrajectory_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Apr  1 08:23:25 EDT 2016
//
// This is a base class for all things trajectory in 3D as a
// function of time.
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"


class TParticleTrajectory
{
  public:
    virtual TVector3D GetPosition (double const&) const = 0;
    virtual TVector3D GetVelocity (double const&) const = 0;
    virtual TVector3D GetAcceleration (double const&) const = 0;

};













#endif

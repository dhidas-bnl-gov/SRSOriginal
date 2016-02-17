#ifndef GUARD_TParticleTrajectory_h
#define GUARD_TParticleTrajectory_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon Dec 14 17:36:45 EST 2015
//
////////////////////////////////////////////////////////////////////


#include <vector>
#include <array>

class TParticleTrajectory
{
  public:
    TParticleTrajectory ();
    ~TParticleTrajectory ();

    void Add (double const, double const, double const, double const, double const, double const);

  private:
    double fZFirstPoint;
    double fZLastPoint;
    int    fZNPoints;
    int    fZNPointsPerMeter;
    double fZStepSize;

    std::vector< std::array<double, 6> > fPositionVelocity;

};




































#endif

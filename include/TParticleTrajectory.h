#ifndef GUARD_TParticleTrajectory_h
#define GUARD_TParticleTrajectory_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon Dec 14 17:36:45 EST 2015
//
////////////////////////////////////////////////////////////////////



class TParticleTrajectory
{
  public:
    TParticleTrajectory ();
    ~TParticleTrajectory ();


  private:
    double fZFirstPoint;
    double fZLastPoint;
    int    fZNPoints;
    int    fZNPointsPerMeter;
    double fZStepSize;

};




































#endif

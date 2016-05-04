#ifndef GUARD_TParticleTrajectory_Points_h
#define GUARD_TParticleTrajectory_Points_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Mar 25 17:18:46 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TParticleTrajectory.h"
#include "TVector3D.h"

#include <vector>


class TParticleTrajectory_Points : public TParticleTrajectory
{
  public:
    TParticleTrajectory_Points (double const&, double const&);
    ~TParticleTrajectory_Points ();

    void AddPoint (double const&, double const&, double const&, double const&, double const&, double const&, double const&, double const&, double const&);

    TVector3D GetPosition (double const&) const;
    TVector3D GetVelocity (double const&) const;
    TVector3D GetAcceleration (double const&) const;

    double GetStartTime () const;

  private:
    double fStartTime;
    double fDeltaTime;

    std::vector<TVector3D*> fX;
    std::vector<TVector3D*> fV;
    std::vector<TVector3D*> fA;

    size_t GetFirstBinForTime (double const&) const;

};




















#endif

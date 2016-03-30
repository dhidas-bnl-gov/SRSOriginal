#ifndef GUARD_TParticleTrajectory_h
#define GUARD_TParticleTrajectory_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Mar 25 17:18:46 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"

#include <vector>


class TParticleTrajectory
{
  public:
    TParticleTrajectory ();
    TParticleTrajectory (double const&, double const&);
    ~TParticleTrajectory ();

    void AddPoint(double const&, double const&, double const&, double const&, double const&, double const&, double const&, double const&, double const&);

    TVector3D GetXAtTime(double const&) const;
    TVector3D GetVAtTime(double const&) const;
    TVector3D GetAAtTime(double const&) const;

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

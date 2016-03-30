#ifndef GUARD_TParticleTrajectory2_h
#define GUARD_TParticleTrajectory2_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Mar 25 17:18:46 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"

#include <vector>


class TParticleTrajectory2
{
  public:
    TParticleTrajectory2 ();
    TParticleTrajectory2 (double const&, double const&);
    ~TParticleTrajectory2 ();

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

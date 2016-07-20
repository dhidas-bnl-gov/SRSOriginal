#ifndef GUARD_TTwiss_h
#define GUARD_TTwiss_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jun  9 15:54:43 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"

class TTwiss
{
  public:
    TTwiss ();
    TTwiss (double const, double const, double const, double const, TVector3D const& X0 = TVector3D(0, 0, 0));
    ~TTwiss ();

    void SetParameters ();
    void SetX0 (TVector3D const&);

  private:

    TVector3D fX0;  // Location that Twiss parameters are defined

    TVector3D fAlpha;
    TVector3D fBeta;
    TVector3D fEpsilon;
    TVector3D fGamma;
};




#endif

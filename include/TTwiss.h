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
    ~TTwiss ();

  private:
    double fTestVariable;

    TVector3D fX0;  // Location that Twiss parameters are defined
};




#endif

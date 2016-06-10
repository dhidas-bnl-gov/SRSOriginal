#ifndef GUARD_TBFieldContainer
#define GUARD_TBFieldContainer
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue May 24 08:59:22 EDT 2016
//
// This class is meant to be a container class for all B-fields in
// a given simulation.  It will sum all B contributions and return
// the sum Bx, By, Bz, or zero where there is no defined field
//
////////////////////////////////////////////////////////////////////

#include <vector>

#include "TBField.h"
#include "TVector3D.h"

class TBFieldContainer
{
  public:
    TBFieldContainer ();
    TBFieldContainer (TBField*);
    ~TBFieldContainer ();

    void AddField (TBField*);

    double    GetBx (double const, double const, double const) const;
    double    GetBy (double const, double const, double const) const;
    double    GetBz (double const, double const, double const) const;
    TVector3D GetB  (double const, double const, double const) const;

  private:
    std::vector<TBField*> fBFields;
};




























#endif
